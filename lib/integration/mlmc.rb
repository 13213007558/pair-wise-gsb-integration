module Integration
  # Result of a multilevel Monte Carlo estimation, see Integration.mlmc.
  #
  # +standard_error+ and +confidence_interval+ quantify the *sampling*
  # (Monte Carlo) error only. They are computed from the sample variances
  # of the per-level difference estimators and say nothing about the
  # discretization bias of the finest level, i.e. the confidence interval
  # is NOT a guarantee that the exact (continuum) expectation lies inside
  # it. +bias_estimate+ is a separate heuristic indicator of the remaining
  # discretization bias and should be inspected alongside the interval.
  class MLMCResult
    # Estimated expectation: sum of the per-level means.
    attr_reader :estimate

    # Standard error of +estimate+ due to Monte Carlo sampling only.
    attr_reader :standard_error

    # Hash mapping level => number of samples actually drawn at that level
    # (pilot plus continuation samples).
    attr_reader :samples

    # Hash mapping level => sample mean of the level difference
    # (level 0 is the plain coarse estimator).
    attr_reader :level_means

    # Hash mapping level => sample variance of the level difference.
    attr_reader :level_variances

    # Hash mapping level => cost per sample used for the budget allocation.
    attr_reader :level_costs

    # Heuristic estimate of the remaining discretization bias
    # |E[P] - E[P_L]|, extrapolated from the two finest level means.
    # This is only an indicator, not a bound.
    attr_reader :bias_estimate

    # Confidence level used for #confidence_interval (default 0.95).
    attr_reader :confidence

    def initialize(attrs = {})
      @estimate        = attrs[:estimate]
      @standard_error  = attrs[:standard_error]
      @samples         = attrs[:samples]
      @level_means     = attrs[:level_means]
      @level_variances = attrs[:level_variances]
      @level_costs     = attrs[:level_costs]
      @bias_estimate   = attrs[:bias_estimate]
      @confidence      = attrs[:confidence] || 0.95
    end

    # Sampling variance of the estimate (square of #standard_error).
    def variance
      standard_error**2
    end

    # Number of levels used (levels are 0..levels-1).
    def levels
      samples.size
    end

    # Total number of samples drawn across all levels.
    def total_samples
      samples.values.inject(0, :+)
    end

    # Normal-approximation confidence interval for the estimate, based
    # solely on the sampling error. It does NOT account for discretization
    # bias; see the class documentation.
    def confidence_interval(level = confidence)
      z = Integration.send(:mlmc_normal_quantile, 0.5 + level / 2.0)
      half = z * standard_error
      [estimate - half, estimate + half]
    end

    def to_h
      {
        estimate: estimate,
        standard_error: standard_error,
        samples: samples.dup,
        level_means: level_means.dup,
        level_variances: level_variances.dup,
        level_costs: level_costs.dup,
        bias_estimate: bias_estimate,
        confidence: confidence
      }
    end
  end

  class << self
    # Multilevel Monte Carlo (MLMC) estimation of an expectation from a
    # hierarchy of discretizations of the same simulation (cheap coarse
    # grids up to expensive fine grids).
    #
    # The caller supplies a *coupled* sampler as a block. The block is
    # called as +sampler.call(level, rng)+ and must return:
    #
    # * for level 0: the coarse estimator P_0 evaluated at one random
    #   input (a numeric, or a one-element array);
    # * for level l >= 1: a two-element array [P_l_fine, P_{l-1}_coarse]
    #   where both values are evaluated at the SAME random input, so that
    #   the difference P_l - P_{l-1} has small variance. Use the passed
    #   +rng+ for all randomness; drawing independent inputs for the fine
    #   and coarse parts destroys the variance reduction.
    #
    # The estimator is
    #
    #   E[P_L] = E[P_0] + sum_{l=1..L} E[P_l - P_{l-1}]
    #
    # with each term averaged over an independently sized sample.
    #
    # Procedure:
    #
    # 1. A pilot of +pilot_samples+ samples is drawn at every level to
    #    estimate the variance V_l of each level difference.
    # 2. The remaining budget is allocated as N_l proportional to
    #    sqrt(V_l / C_l), which minimizes the total sampling variance for
    #    a fixed total cost (or minimizes cost for a fixed target
    #    variance). Only the shortfall over the pilot samples is drawn.
    # 3. Levels are appended while the extrapolated remaining bias
    #    |mean_L| / (|mean_{L-1}| / |mean_L| - 1) exceeds the target
    #    sampling error (i.e. the finest-level correction has not yet
    #    decayed enough for the discretization bias to be negligible at
    #    the requested accuracy), up to +max_levels+. When the ratio of
    #    consecutive level means does not indicate decay, |mean_L| itself
    #    is used as the bias indicator.
    #
    # Options:
    # [:tolerance]     Target for the sampling standard error of the
    #                  estimate. Mutually exclusive with +:budget+.
    # [:budget]        Total cost budget (sum over levels of
    #                  cost_per_sample * samples). Mutually exclusive with
    #                  +:tolerance+.
    # [:costs]         Cost per sample at each level: a lambda taking the
    #                  level, or a Hash/Array indexed by level.
    #                  Default: ->(l) { 2.0**l }.
    # [:pilot_samples] Samples per level in the pilot round. Default: 25.
    # [:min_levels]    Number of levels (0..min_levels-1) used initially.
    #                  Default: 3.
    # [:max_levels]    Maximum number of levels. Default: 10.
    # [:rng]           Random source passed to the sampler.
    #                  Default: Random.new.
    # [:confidence]    Confidence level for MLMCResult#confidence_interval.
    #                  Default: 0.95.
    #
    # Returns an MLMCResult with the expectation estimate, the sampling
    # standard error and the per-level sample counts.
    #
    # NOTE: the returned standard error and confidence interval reflect
    # Monte Carlo sampling error only. They do not include the
    # discretization bias of the finest level; check
    # MLMCResult#bias_estimate separately.
    def mlmc(options = {}, &sampler)
      fail 'No sampler passed' unless block_given?

      tolerance = options[:tolerance]
      budget    = options[:budget]
      if tolerance.nil? && budget.nil?
        fail 'Pass either :tolerance (target sampling error) or :budget (total cost)'
      end
      fail 'Pass either :tolerance or :budget, not both' if tolerance && budget

      pilot      = options.fetch(:pilot_samples, 25).to_i
      fail ':pilot_samples must be >= 2' if pilot < 2
      min_levels = options.fetch(:min_levels, 3).to_i
      max_levels = options.fetch(:max_levels, 10).to_i
      fail ':min_levels must be >= 2' if min_levels < 2
      fail ':max_levels must be >= :min_levels' if max_levels < min_levels
      rng        = options[:rng] || Random.new
      confidence = options.fetch(:confidence, 0.95)
      cost_fn    = mlmc_cost_function(options[:costs])

      stats = {} # level => {n:, mean:, m2:}
      n_levels = min_levels

      loop do
        (0...n_levels).each do |l|
          mlmc_top_up(stats, l, pilot, rng, sampler)
        end

        allocation = mlmc_allocate(stats, n_levels, cost_fn, tolerance, budget)
        target_se  = tolerance || Math.sqrt(mlmc_total_variance(stats, allocation))

        if n_levels < max_levels &&
           mlmc_bias_indicator(stats, n_levels) > target_se
          n_levels += 1
          next
        end

        (0...n_levels).each do |l|
          mlmc_top_up(stats, l, allocation[l], rng, sampler)
        end
        break
      end

      estimate = 0.0
      means = {}
      variances = {}
      costs = {}
      samples = {}
      (0...n_levels).each do |l|
        s = stats[l]
        means[l]     = s[:mean]
        variances[l] = s[:n] > 1 ? s[:m2] / (s[:n] - 1.0) : 0.0
        costs[l]     = cost_fn.call(l)
        samples[l]   = s[:n]
        estimate    += s[:mean]
      end

      MLMCResult.new(
        estimate: estimate,
        standard_error: Math.sqrt(mlmc_total_variance(stats, samples)),
        samples: samples,
        level_means: means,
        level_variances: variances,
        level_costs: costs,
        bias_estimate: mlmc_bias_indicator(stats, n_levels),
        confidence: confidence
      )
    end

    private

    # Draw samples at +level+ until stats[level][:n] reaches +target+.
    def mlmc_top_up(stats, level, target, rng, sampler)
      s = stats[level] ||= { n: 0, mean: 0.0, m2: 0.0 }
      while s[:n] < target
        y = mlmc_draw(level, rng, sampler)
        s[:n] += 1
        delta = y - s[:mean]
        s[:mean] += delta / s[:n]
        s[:m2] += delta * (y - s[:mean])
      end
    end

    # One coupled sample of the level difference (or of P_0 at level 0).
    def mlmc_draw(level, rng, sampler)
      out = sampler.call(level, rng)
      if level.zero?
        v = out.is_a?(Array) ? out[0] : out
        fail "Sampler returned a non-numeric value at level 0: #{out.inspect}" unless v.is_a?(Numeric)
        v.to_f
      else
        unless out.is_a?(Array) && out.size == 2 && out.all? { |o| o.is_a?(Numeric) }
          fail "Sampler at level #{level} must return [fine_value, coarse_value] " \
               "computed from the same random input, got: #{out.inspect}"
        end
        out[0].to_f - out[1].to_f
      end
    end

    def mlmc_cost_function(costs)
      case costs
      when nil
        ->(l) { 2.0**l }
      when Proc
        costs
      when Hash, Array
        ->(l) { costs[l] }
      else
        fail ':costs must be a lambda, a Hash or an Array'
      end
    end

    # Optimal MLMC sample allocation N_l proportional to sqrt(V_l / C_l),
    # scaled to meet either a target sampling standard error or a total
    # cost budget. Never below the number of samples already drawn.
    def mlmc_allocate(stats, n_levels, cost_fn, tolerance, budget)
      v = {}
      c = {}
      sum_vc = 0.0
      (0...n_levels).each do |l|
        s = stats[l]
        v[l] = s[:n] > 1 ? s[:m2] / (s[:n] - 1.0) : 0.0
        c[l] = cost_fn.call(l).to_f
        fail "Cost at level #{l} must be positive" unless c[l] > 0
        sum_vc += Math.sqrt(v[l] * c[l])
      end

      allocation = {}
      if sum_vc.zero?
        (0...n_levels).each { |l| allocation[l] = stats[l][:n] }
        return allocation
      end

      (0...n_levels).each do |l|
        n = if tolerance
              sum_vc * Math.sqrt(v[l] / c[l]) / (tolerance**2)
            else
              budget * Math.sqrt(v[l] / c[l]) / sum_vc
            end
        allocation[l] = [n.ceil, stats[l][:n]].max
      end
      allocation
    end

    def mlmc_total_variance(stats, allocation)
      total = 0.0
      allocation.each do |l, n|
        s = stats[l]
        next if n.nil? || n < 1
        var = s[:n] > 1 ? s[:m2] / (s[:n] - 1.0) : 0.0
        total += var / n
      end
      total
    end

    # Heuristic indicator of the remaining discretization bias after the
    # finest level, extrapolated from the two finest level means assuming
    # geometric decay of the corrections. Falls back to the finest level
    # mean itself when no decay is observable. Not a rigorous bound.
    def mlmc_bias_indicator(stats, n_levels)
      finest = stats[n_levels - 1][:mean].abs
      return finest if n_levels < 2

      previous = stats[n_levels - 2][:mean].abs
      return finest if finest.zero? || previous <= finest

      ratio = previous / finest
      finest / (ratio - 1.0)
    end

    # Inverse of the standard normal CDF (Acklam's rational approximation,
    # max absolute error about 1.15e-9 in the central region).
    def mlmc_normal_quantile(p)
      fail 'p must be in (0, 1)' unless p > 0.0 && p < 1.0

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

      plow = 0.02425
      phigh = 1 - plow

      if p < plow
        q = Math.sqrt(-2 * Math.log(p))
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
          ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
      elsif p > phigh
        q = Math.sqrt(-2 * Math.log(1 - p))
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
          ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
      else
        q = p - 0.5
        r = q * q
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
          (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
      end
    end
  end
end
