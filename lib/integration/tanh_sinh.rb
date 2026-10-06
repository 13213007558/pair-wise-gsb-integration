module Integration
  # Result of a tanh-sinh (double exponential) quadrature, as returned by
  # Integration.tanh_sinh.
  class TanhSinhResult
    # Estimated value of the integral.
    attr_reader :value

    # Discretization (mesh) error estimate: absolute difference between the
    # two finest mesh levels.
    attr_reader :error

    # Transformation tail truncation estimate: magnitude of the first
    # neglected terms at both ends of the transformed domain.
    attr_reader :tail_error

    # Number of integrand evaluations performed.
    attr_reader :evaluations

    # Number of mesh refinements performed (0 means only the initial mesh
    # was evaluated).
    attr_reader :levels

    # Why the iteration stopped. One of:
    # [:ok]          both the discretization error and the tail truncation
    #                estimate are within the requested tolerance.
    # [:max_levels]  the refinement limit was reached before the tolerance
    #                was met; the reported value was NOT verified against
    #                the requested tolerance.
    # [:divergent]   error estimates grow under mesh refinement (or the
    #                partial sums overflow); the integral appears to be
    #                divergent.
    # [:non_finite]  the integrand evaluated to NaN at an interior point.
    attr_reader :status

    def initialize(value, error, tail_error, evaluations, levels, status)
      @value = value
      @error = error
      @tail_error = tail_error
      @evaluations = evaluations
      @levels = levels
      @status = status
    end

    # True only when the requested tolerance was actually met (+:ok+).
    def converged?
      @status == :ok
    end
  end

  class << self
    # Default options for tanh-sinh quadrature.
    TANH_SINH_DEFAULTS = {
      tolerance: 1e-10,   # absolute tolerance for both error estimates
      initial_step: 0.5,  # initial mesh size h in the transformed coordinate
      max_levels: 10,     # maximum number of mesh halvings
      max_points: 1 << 16 # safety cap on points per side per level
    }.freeze

    # Tanh-sinh (double exponential) quadrature over a finite interval.
    #
    # Well suited for integrands with integrable endpoint singularities
    # such as +log(x)+ or +1/sqrt(x)+: the integrand is never evaluated at
    # the endpoints, so the caller does not need to provide finite (or any)
    # function values there.
    #
    # The mesh (discretization) error and the transformation tail
    # truncation are estimated independently; convergence is only reported
    # when BOTH are within +tolerance+. Agreement between two successive
    # mesh levels alone is never taken as proof of convergence.
    #
    # Options are:
    # [:tolerance]    Absolute tolerance applied to both the discretization
    #                 error and the tail truncation estimate.
    #                 Default: 1e-10.
    # [:initial_step] Initial mesh size in the transformed coordinate.
    #                 Default: 0.5.
    # [:max_levels]   Maximum number of mesh halvings. Default: 10.
    # [:max_points]   Safety cap on quadrature points per side per level.
    #
    # Returns a TanhSinhResult; check #converged? / #status before trusting
    # #value. A genuinely divergent integral is reported with a non-+:ok+
    # status, never as tolerance achieved.
    def tanh_sinh(t1, t2, options = {}, &f)
      fail 'No function passed' unless block_given?
      fail 'Non-numeric bounds' unless t1.is_a?(Numeric) && t2.is_a?(Numeric)
      fail 'tanh_sinh does not support infinite bounds' if infinite?(t1) || infinite?(t2)
      fail 'Zero-length integration interval' if t1 == t2

      opts = TANH_SINH_DEFAULTS.merge(options)
      tolerance = opts[:tolerance].to_f
      max_levels = opts[:max_levels].to_i
      max_points = opts[:max_points].to_i
      fail 'Tolerance must be positive' unless tolerance > 0.0
      fail 'initial_step must be positive' unless opts[:initial_step].to_f > 0.0

      sign = 1.0
      lower = t1.to_f
      upper = t2.to_f
      if lower > upper
        lower, upper = upper, lower
        sign = -1.0
      end
      h = opts[:initial_step].to_f
      evaluations = 0

      sum, tail_error, ev, failure = tanh_sinh_level(lower, upper, h, 0, 1,
                                                     tolerance, max_points, f)
      evaluations += ev
      if failure
        return TanhSinhResult.new(sign * sum, Float::INFINITY, tail_error,
                                  evaluations, 0, failure)
      end

      levels = 0
      disc_error = Float::INFINITY
      consecutive_increases = 0
      slow_tail_levels = 0
      previous_tail_error = tail_error
      status = nil

      loop do
        h *= 0.5
        odd_sum, tail_error, ev, failure = tanh_sinh_level(lower, upper, h, 1, 2,
                                                           tolerance, max_points, f)
        evaluations += ev
        if failure
          sum = odd_sum
          status = failure
          break
        end

        refined = 0.5 * sum + odd_sum
        previous_disc_error = disc_error
        disc_error = (refined - sum).abs
        sum = refined
        levels += 1

        if !sum.finite? || !disc_error.finite?
          status = :divergent
          break
        elsif disc_error <= tolerance && tail_error <= tolerance
          status = :ok
          break
        elsif levels >= max_levels
          status = :max_levels
          break
        end

        # Divergence check: the discretization error must shrink under
        # refinement. Two consecutive increases above the tolerance mean the
        # partial sums are blowing up (e.g. a non-integrable singularity).
        if disc_error > previous_disc_error
          consecutive_increases += 1
          if consecutive_increases >= 2
            status = :divergent
            break
          end
        else
          consecutive_increases = 0
        end

        # Divergence check on the transformation tail: for an integrable
        # function the neglected tail decays double exponentially once the
        # truncation point enters the asymptotic region. A tail that stays
        # above the tolerance while decaying only geometrically (less than
        # a factor of 4 per level) over several consecutive levels means
        # the transformed integrand itself does not decay: the integral is
        # divergent (e.g. a non-integrable singularity like 1/x).
        if tail_error > tolerance && tail_error > 0.25 * previous_tail_error
          slow_tail_levels += 1
          if slow_tail_levels >= 4
            status = :divergent
            break
          end
        else
          slow_tail_levels = 0
        end
        previous_tail_error = tail_error
      end

      disc_error = Float::INFINITY if disc_error.nan?
      TanhSinhResult.new(sign * sum, disc_error, tail_error,
                         evaluations, levels, status)
    end

    private

    # Transformed abscissa/weight pair for tanh-sinh quadrature.
    #
    #   x(t) = center + half * tanh(pi/2 * sinh(t))
    #   w(t) = half * (pi/2) * cosh(t) / cosh(pi/2 * sinh(t))^2
    #
    # The abscissa is computed as a distance from the nearest endpoint via
    # 1 / (1 + exp(2*|arg|)); the naive center + half*tanh(arg) form loses
    # all precision to catastrophic cancellation exactly near the endpoints,
    # where the integrable singularities this method targets live.
    #
    # Returns a weight of exactly 0.0 whenever the abscissa coincides with
    # an endpoint or the weight underflows, so the caller never evaluates
    # the integrand at the endpoints.
    def tanh_sinh_point(lower, upper, t)
      half_pi = Math::PI / 2
      arg = half_pi * Math.sinh(t)
      return [arg > 0 ? upper : lower, 0.0] if arg.abs > 350.0

      width = upper - lower
      if arg >= 0.0
        x = upper - width / (1.0 + Math.exp(2.0 * arg))
      else
        x = lower + width / (1.0 + Math.exp(-2.0 * arg))
      end
      return [x, 0.0] if x <= lower || x >= upper

      cosh_arg = Math.cosh(arg)
      weight = 0.5 * width * half_pi * Math.cosh(t) / (cosh_arg * cosh_arg)
      weight = 0.0 unless weight.finite?
      [x, weight]
    end

    # Trapezoidal sum (times +h+) over the k-indices +k_start, k_start +
    # k_step, ...+ on both sides of 0. For +k_start == 0+ the center point
    # is included once and both sides start at +k_step+.
    #
    # Each side is truncated independently once a term's contribution falls
    # below a quarter of the tolerance or the weight underflows; the
    # magnitude of the first neglected term is reported as the tail
    # truncation estimate (double exponential decay makes the remaining
    # tail comparable to it).
    #
    # Returns [sum, tail_estimate, evaluations, failure_status].
    def tanh_sinh_level(lower, upper, h, k_start, k_step, tolerance, max_points, f)
      evaluations = 0
      sum = 0.0
      tail = 0.0

      if k_start == 0
        x, weight = tanh_sinh_point(lower, upper, 0.0)
        term = weight * f.call(x)
        evaluations += 1
        return [sum, tail, evaluations, :non_finite] if term.nan?
        sum += term
        k_start = k_step
      end

      [1, -1].each do |side|
        k = k_start
        last_scaled = 0.0
        while k <= max_points
          x, weight = tanh_sinh_point(lower, upper, side * k * h)
          break if weight <= 0.0 # weight underflowed: remaining tail negligible

          term = weight * f.call(x)
          evaluations += 1
          return [h * sum, tail, evaluations, :non_finite] if term.nan?

          last_scaled = h * term.abs
          break if last_scaled <= 0.25 * tolerance

          sum += term
          k += k_step
        end
        tail += last_scaled
      end

      [h * sum, tail, evaluations, nil]
    end
  end
end
