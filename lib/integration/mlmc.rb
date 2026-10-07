module Integration
  # Result object returned by Integration.mlmc.
  #
  # Attributes:
  # [estimate]        Multilevel Monte Carlo estimate of the expectation.
  # [standard_error]  Standard error of the estimate, i.e.
  #                   sqrt(sum_l Var(Y_l) / N_l). This reflects Monte Carlo
  #                   sampling error ONLY; see #confidence_interval.
  # [samples]         Hash mapping each level to the number of coupled
  #                   samples actually drawn on that level (pilot included).
  # [level_means]     Hash mapping each level to the sample mean of its
  #                   level difference Y_l.
  # [level_variances] Hash mapping each level to the (unbiased) sample
  #                   variance of its level difference Y_l.
  # [bias_estimate]   Heuristic estimate of the remaining discretization
  #                   bias E[P] - E[P_L], extrapolated from the decay of the
  #                   level-difference means. It is only as reliable as the
  #                   geometric-decay assumption behind it.
  # [cost_spent]      Total cost units consumed (sum of samples * cost).
  class MLMCResult
    attr_reader :estimate, :standard_error, :samples, :level_means,
                :level_variances, :bias_estimate, :cost_spent

    def initialize(attributes)
      @estimate        = attributes[:estimate]
      @standard_error  = attributes[:standard_error]
      @samples         = attributes[:samples]
      @level_means     = attributes[:level_means]
      @level_variances = attributes[:level_variances]
      @bias_estimate   = attributes[:bias_estimate]
      @cost_spent      = attributes[:cost_spent]
    end

    # Levels that were sampled, in ascending order.
    def levels
      @samples.keys.sort
    end

    # Normal-approximation confidence interval for the estimate.
    #
    # WARNING: this interval quantifies Monte Carlo *sampling* error only.
    # It is NOT a strict guarantee that the exact expectation of the
    # underlying (continuous) model lies inside it: the finest level still
    # carries a discretization (bias) error, which #bias_estimate merely
    # approximates under a heuristic geometric-decay assumption. Treat the
    # interval as "estimate +/- sampling noise", not as a bound on the total
    # error.
    def confidence_interval(confidence = 0.95)
      fail 'confidence must be in (0, 1)' unless confidence > 0 && confidence < 1
      z = Integration.send(:normal_quantile, 0.5 + confidence / 2.0)
      half_width = z * @standard_error
      [@estimate - half_width, @estimate + half_width]
    end
  end

  class << self
    # Multilevel Monte Carlo (MLMC) estimation of an expectation E[P], where
    # P is approximated by a hierarchy of discretizations P_0, P_1, ...
    # (e.g. the same simulation on increasingly fine grids, cheap coarse
    # versions and expensive fine ones).
    #
    # The telescoping identity E[P_L] = E[P_0] + sum_{l=1..L} E[P_l - P_{l-1}]
    # is used: each level difference Y_l = P_l - P_{l-1} is estimated by
    # Monte Carlo with N_l samples, and most samples are spent on the cheap
    # coarse levels where the differences have the largest variance.
    #
    # Options:
    # [:sampler] REQUIRED. Callable +sampler.call(level, n, rng)+ returning an
    #            Array of +n+ i.i.d. samples of the level difference
    #            Y_l = P_l(omega) - P_{l-1}(omega) (Y_0 = P_0(omega)).
    #            Crucially, the coarse and fine approximations inside one
    #            sample must share the SAME random input +omega+ (coupling);
    #            otherwise the level differences do not cancel and the
    #            method degenerates to plain Monte Carlo. +rng+ is the
    #            Random instance passed via the +:rng+ option, so callers
    #            keep full control of the random source.
    # [:costs]   REQUIRED. Cost of one coupled sample per level, as a
    #            callable +costs.call(level)+ or a Hash/Array indexed by
    #            level. Used for optimal budget allocation.
    # [:budget]  REQUIRED. Total cost units available. The pilot phase
    #            always runs (even if it alone exceeds the budget); the
    #            remaining budget is then distributed over the levels.
    # [:pilot_samples] Samples per level in the pilot phase. Default: 25.
    #            Must be >= 2 so level variances are estimable.
    # [:levels]  Initial set of levels. Default: [0, 1].
    # [:extend_levels] Whether finer levels may be appended after the pilot
    #            phase. Default: true. See the note on level extension below.
    # [:max_level] Hard cap on the finest level. Default: 10.
    # [:tolerance] Optional target for the total root-mean-square error.
    #            When given, level extension stops once the bias proxy drops
    #            below tolerance/sqrt(2). When omitted, the current sampling
    #            standard error is used as the allowed bias instead, so that
    #            levels are added while the (estimated) discretization bias
    #            would dominate the sampling error.
    # [:rng]     Random source handed to the sampler. Default: Random.new.
    #            Pass Random.new(seed) for reproducible runs.
    #
    # Pilot phase and budget allocation:
    # A small pilot batch of +:pilot_samples+ samples is drawn on every
    # current level to estimate the level-difference variances V_l. The
    # remaining budget B is then distributed to minimize the estimator
    # variance sum_l V_l / N_l subject to sum_l C_l N_l <= B, giving the
    # standard optimal allocation N_l proportional to sqrt(V_l / C_l).
    # Leftover cost from integer rounding is spent greedily on the level
    # with the largest variance reduction per unit cost.
    #
    # Level extension rationale:
    # Assuming the level-difference means decay geometrically,
    # |E[Y_l]| ~ 2^(-alpha*l), the bias remaining after the finest level L
    # is approximately |E[Y_L]| / (2^alpha - 1), with alpha estimated from
    # the ratio of the last two level means (defaulting to 1 when the ratio
    # is unusable). A finer level is appended while this proxy exceeds the
    # allowed bias (see +:tolerance+), while +:max_level+ is not reached,
    # and while the budget still affords a pilot batch on the new level.
    # This is a heuristic based on the observed decay of the level means,
    # not a proof about the true discretization error.
    #
    # Returns an Integration::MLMCResult. Note that the reported standard
    # error and confidence interval cover sampling error only; they do not
    # strictly bound the discretization bias of the finest level.
    def mlmc(options = {})
      fail 'Missing coupled sampler (:sampler)' unless options[:sampler]
      fail 'Missing level costs (:costs)' unless options[:costs]
      fail 'Missing total budget (:budget)' unless options[:budget]

      sampler = options[:sampler]
      costs   = options[:costs]
      cost_of = costs.respond_to?(:call) ? costs : ->(level) { costs[level] }
      budget  = options[:budget].to_f
      pilot   = (options[:pilot_samples] || 25).to_i
      fail ':pilot_samples must be >= 2' if pilot < 2
      levels    = (options[:levels] || [0, 1]).sort
      extend    = options.fetch(:extend_levels, true)
      max_level = options[:max_level] || 10
      tolerance = options[:tolerance]
      rng       = options[:rng] || Random.new

      values    = {}
      means     = {}
      variances = {}
      samples   = Hash.new(0)
      spent     = 0.0

      draw = lambda do |level, n|
        draws = sampler.call(level, n, rng)
        (values[level] ||= []).concat(draws)
        samples[level] += draws.size
        spent += draws.size * cost_of.call(level)
      end

      recompute = lambda do
        levels.each do |level|
          v = values[level]
          mean = v.inject(0.0, :+) / v.size
          means[level] = mean
          variances[level] = v.inject(0.0) { |acc, x| acc + (x - mean)**2 } / (v.size - 1)
        end
      end

      # Pilot phase: small batch on every level to estimate V_l.
      levels.each { |level| draw.call(level, pilot) }

      bias_estimate = nil
      if extend
        loop do
          recompute.call
          finest = levels.max
          break if finest >= max_level
          alpha = mlmc_weak_rate(means, levels)
          bias_estimate = means[finest].abs / ((2.0**alpha) - 1.0)
          standard_error = Math.sqrt(levels.inject(0.0) { |acc, l| acc + variances[l] / samples[l] })
          allowed_bias = tolerance ? tolerance / Math.sqrt(2.0) : standard_error
          break unless bias_estimate > allowed_bias
          break if spent + pilot * cost_of.call(finest + 1) > budget
          levels << finest + 1
          draw.call(finest + 1, pilot)
        end
      end
      recompute.call

      # Budget allocation: N_l proportional to sqrt(V_l / C_l).
      remaining = budget - spent
      if remaining > 0
        mlmc_allocate(levels, variances, samples, cost_of, remaining).each do |level, n|
          draw.call(level, n) if n > 0
        end
        recompute.call
      end

      estimate = levels.inject(0.0) { |acc, l| acc + means[l] }
      sampling_variance = levels.inject(0.0) { |acc, l| acc + variances[l] / samples[l] }
      if bias_estimate.nil?
        alpha = mlmc_weak_rate(means, levels)
        bias_estimate = means[levels.max].abs / ((2.0**alpha) - 1.0)
      end

      MLMCResult.new(
        estimate: estimate,
        standard_error: Math.sqrt(sampling_variance),
        samples: levels.each_with_object({}) { |l, h| h[l] = samples[l] },
        level_means: means.select { |l, _| levels.include?(l) },
        level_variances: variances.select { |l, _| levels.include?(l) },
        bias_estimate: bias_estimate,
        cost_spent: spent
      )
    end

    private

    # Estimate the weak convergence rate alpha of the level-difference
    # means, assuming |E[Y_l]| ~ 2^(-alpha*l). Falls back to alpha = 1 when
    # the observed means are zero or do not decay.
    def mlmc_weak_rate(means, levels)
      sorted = levels.sort
      fine = sorted[-1]
      coarse = sorted[-2]
      return 1.0 if coarse.nil?
      m_fine = means[fine].abs
      m_coarse = means[coarse].abs
      return 1.0 if m_fine.zero? || m_coarse.zero?
      rate = (m_coarse / m_fine)**(1.0 / (fine - coarse))
      rate > 1.0 ? Math.log2(rate) : 1.0
    end

    # Distribute +budget_left+ cost units over the levels. The bulk follows
    # the variance-optimal allocation N_l proportional to sqrt(V_l / C_l);
    # the integer-rounding leftover is spent greedily on the sample with the
    # largest variance reduction per unit cost.
    def mlmc_allocate(levels, variances, samples, cost_of, budget_left)
      extra = Hash.new(0)
      weights = {}
      levels.each do |level|
        variance = variances[level]
        cost = cost_of.call(level)
        weights[level] = (variance > 0 && cost > 0) ? Math.sqrt(variance / cost) : 0.0
      end
      denom = levels.inject(0.0) { |acc, l| acc + cost_of.call(l) * weights[l] }
      if denom > 0
        levels.each { |l| extra[l] = (budget_left * weights[l] / denom).floor }
      end

      leftover = budget_left - levels.inject(0.0) { |acc, l| acc + extra[l] * cost_of.call(l) }
      loop do
        best = nil
        best_gain = 0.0
        levels.each do |level|
          cost = cost_of.call(level)
          next if cost <= 0 || cost > leftover
          n = samples[level] + extra[level]
          gain = variances[level] / (cost * n * (n + 1.0))
          if gain > best_gain
            best_gain = gain
            best = level
          end
        end
        break if best.nil?
        extra[best] += 1
        leftover -= cost_of.call(best)
      end
      extra
    end

    # Inverse of the standard normal CDF (Acklam's rational approximation).
    def normal_quantile(p)
      a = [-3.969683028665376e+01, 2.209460984245205e+02,
           -2.759285104469687e+02, 1.383577518672690e+02,
           -3.066479806614716e+01, 2.506628277459239e+00]
      b = [-5.447609879822406e+01, 1.615858368580409e+02,
           -1.556989798598866e+02, 6.680131188771972e+01,
           -1.328068155288572e+01]
      c = [-7.784894002430293e-03, -3.223964580411365e-01,
           -2.400758277161838e+00, -2.549732539343734e+00,
           4.374664141464968e+00, 2.938163982698783e+00]
      d = [7.784695709041462e-03, 3.224671290700398e-01,
           2.445134137142996e+00, 3.754408661907416e+00]
      p_low = 0.02425
      p_high = 1 - p_low

      if p < p_low
        q = Math.sqrt(-2 * Math.log(p))
        num = c[0]; den = d[0]
        1.upto(5) { |i| num = num * q + c[i] }
        1.upto(3) { |i| den = den * q + d[i] }
        num / (den * q + 1)
      elsif p > p_high
        q = Math.sqrt(-2 * Math.log(1 - p))
        num = c[0]; den = d[0]
        1.upto(5) { |i| num = num * q + c[i] }
        1.upto(3) { |i| den = den * q + d[i] }
        -num / (den * q + 1)
      else
        q = p - 0.5
        r = q * q
        num = a[0]; den = b[0]
        1.upto(5) { |i| num = num * r + a[i] }
        1.upto(4) { |i| den = den * r + b[i] }
        num * q / (den * r + 1)
      end
    end
  end
end
