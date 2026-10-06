module Integration
  class ValidatedResult
    attr_reader :lower_bound, :upper_bound, :error, :state,
                :evaluations, :interval_count

    def initialize(lower_bound, upper_bound, state, evaluations, interval_count)
      @lower_bound = lower_bound
      @upper_bound = upper_bound
      @state = state
      @evaluations = evaluations
      @interval_count = interval_count

      @error = if lower_bound.respond_to?(:infinite?) &&
                   (!lower_bound.finite? || !upper_bound.finite?)
        upper_bound - lower_bound
      else
        Rational(upper_bound) - Rational(lower_bound)
      end
    end

    def lower
      lower_bound
    end

    def upper
      upper_bound
    end

    def low
      lower_bound
    end

    def high
      upper_bound
    end

    def width
      error
    end

    def intervals
      interval_count
    end

    def converged?
      state == :converged
    end

    def budget_exhausted?
      state == :budget_exhausted
    end

    def status
      state
    end

    def target_reached?
      converged?
    end

    def reached_target?
      converged?
    end

    def unrefinable_interval?
      state == :unrefinable_interval
    end

    def complete?
      converged?
    end

    def to_a
      [lower_bound, upper_bound]
    end

    def bounds
      to_a
    end

    def deconstruct
      to_a
    end

    def [](index)
      return to_a[index] if index.is_a?(Integer)
      to_h[index] || to_h[index.to_sym]
    end

    def to_h
      {
        lower_bound: lower_bound,
        upper_bound: upper_bound,
        error: error,
        state: state,
        evaluations: evaluations,
        interval_count: interval_count
      }
    end
  end

  class << self
    # Adaptively integrate an interval-valued function oracle.
    #
    # The block must return [lower, upper] (or a Range) such that every value
    # of the integrand on the requested subinterval lies between those bounds.
    # Integral bounds are accumulated with exact Rational arithmetic. Float
    # bounds are accepted only when +outward_round+ is supplied; use +:float+
    # for one-ULP outward rounding or pass a callable taking (value,
    # :down/:up).
    #
    # Returns a ValidatedResult. Its +state+ is +:converged+ when the requested
    # width is reached, otherwise +:budget_exhausted+ or
    # +:unrefinable_interval+; the current validated bounds are returned in
    # every state.
    def validated_integrate(t1, t2, options = {}, &oracle)
      fail ArgumentError, 'An enclosure oracle block is required' unless block_given?
      fail ArgumentError, 'Bounds must be numeric' unless t1.is_a?(Numeric) && t2.is_a?(Numeric)
      fail ArgumentError, 'Infinite bounds are not supported' unless finite_number?(t1) && finite_number?(t2)

      return validated_zero_result(t1, options) if t1 == t2

      orientation = t1 <= t2 ? 1 : -1
      left_input = [t1, t2].min
      right_input = [t1, t2].max
      rounder = validated_rounder(options[:outward_round] || options[:rounding])

      if left_input.is_a?(Float) || right_input.is_a?(Float)
        fail ArgumentError, 'Float bounds require an explicit outward_round mode' if rounder.nil?
      end

      tolerance = exact_number(options.fetch(:tolerance, Rational(1, 10**10)))
      max_evaluations = options[:max_evaluations] || options[:budget] ||
                        options[:evaluation_budget] || options[:max_oracle_calls] ||
                        10_000

      fail ArgumentError, 'max_evaluations must be a positive Integer' unless max_evaluations.is_a?(Integer) && max_evaluations > 0
      fail ArgumentError, 'tolerance must be a finite non-negative number' unless finite_number?(tolerance)
      fail ArgumentError, 'tolerance must be non-negative' unless tolerance >= 0

      leaves = []
      evaluations = 0

      add_leaf = proc do |a_value, b_value|
        low_value, high_value = validated_enclosure(oracle.call(a_value, b_value))
        width = exact_number(b_value) - exact_number(a_value)
        leaves << {
          a_value: a_value,
          b_value: b_value,
          a: exact_number(a_value),
          b: exact_number(b_value),
          low: low_value,
          high: high_value,
          error: width * (high_value - low_value)
        }
      end

      add_leaf.call(left_input, right_input)
      evaluations = 1

      state = :budget_exhausted if evaluations >= max_evaluations

      while state.nil?
        leaf = leaves.max_by { |candidate| candidate[:error] }
        split_value = validated_split_point(leaf, rounder)

        if split_value.nil?
          state = :unrefinable_interval
          break
        end

        if evaluations + 2 > max_evaluations
          state = :budget_exhausted
          break
        end

        leaves.delete(leaf)
        add_leaf.call(leaf[:a_value], split_value)
        add_leaf.call(split_value, leaf[:b_value])
        evaluations += 2

        state = :budget_exhausted if evaluations >= max_evaluations
      end

      build_validated_result(leaves, state, tolerance, rounder, evaluations, orientation)
    end

    private

    def finite_number?(value)
      return false unless value.is_a?(Numeric)
      return false if value.respond_to?(:real?) && !value.real?
      return false if value.respond_to?(:nan?) && value.nan?
      !(value.respond_to?(:infinite?) && value.infinite?)
    end

    def exact_number(value)
      return value if value.is_a?(Rational)
      Rational(value)
    end

    def validated_rounder(spec)
      return nil if spec.nil?
      return method(:float_outward_round) if spec == :float
      return spec if spec.respond_to?(:call)

      fail ArgumentError, 'outward_round must be :float or a callable object'
    end

    def float_outward_round(value, direction)
      nearest = value.to_f

      if nearest.finite?
        return nearest.prev_float if direction == :down

        nearest.next_float
      elsif nearest > 0
        direction == :down ? Float::MAX : Float::INFINITY
      else
        direction == :down ? -Float::INFINITY : -Float::MAX
      end
    end

    def apply_outward_rounding(rounder, value, direction)
      rounded = rounder.call(value, direction)

      fail ArgumentError, 'Outward rounding must return a Numeric value' unless rounded.is_a?(Numeric)
      fail ArgumentError, 'Outward rounding returned NaN' if rounded.respond_to?(:nan?) && rounded.nan?
      if !(rounded.respond_to?(:infinite?) && rounded.infinite?)
        rounded_value = exact_number(rounded)
        valid = direction == :down ? rounded_value <= value : rounded_value >= value
        fail ArgumentError, "Rounding did not round #{direction}" unless valid
      end

      rounded
    end

    def validated_split_point(leaf, rounder)
      midpoint = (leaf[:a] + leaf[:b]) / 2

      if rounder.nil?
        return midpoint
      end

      [:up, :down].each do |direction|
        candidate = apply_outward_rounding(rounder, midpoint, direction)
        next unless finite_number?(candidate)

        candidate_value = exact_number(candidate)
        return candidate if candidate_value > leaf[:a] && candidate_value < leaf[:b]
      end

      nil
    end

    def validated_enclosure(enclosure)
      values = enclosure.is_a?(Range) ? [enclosure.first, enclosure.last] : Array(enclosure)

      unless values.length == 2 && values.all? { |value| finite_number?(value) }
        fail ArgumentError, 'Oracle must return a two-sided numeric enclosure [lower, upper]'
      end

      lower = exact_number(values[0])
      upper = exact_number(values[1])

      unless lower <= upper
        fail ArgumentError, 'Oracle enclosure lower bound must not exceed upper bound'
      end

      [lower, upper]
    end

    def validated_zero_result(bound, _options)
      zero = bound.is_a?(Float) ? 0.0 : Rational(0)
      ValidatedResult.new(zero, zero, :converged, 0, 1)
    end

    def build_validated_result(leaves, current_state, tolerance, rounder, evaluations, orientation)
      exact_lower = leaves.inject(Rational(0)) do |sum, leaf|
        sum + (leaf[:b] - leaf[:a]) * leaf[:low]
      end
      exact_upper = leaves.inject(Rational(0)) do |sum, leaf|
        sum + (leaf[:b] - leaf[:a]) * leaf[:high]
      end

      if orientation < 0
        exact_lower, exact_upper = -exact_upper, -exact_lower
      end

      if rounder.nil?
        lower = exact_lower
        upper = exact_upper
      else
        lower = apply_outward_rounding(rounder, exact_lower, :down)
        upper = apply_outward_rounding(rounder, exact_upper, :up)
      end

      width = if lower.respond_to?(:finite?) && upper.respond_to?(:finite?) && lower.finite? && upper.finite?
        exact_number(upper) - exact_number(lower)
      else
        upper - lower
      end
      state = width <= tolerance ? :converged : current_state

      ValidatedResult.new(lower, upper, state, evaluations, leaves.length)
    end
  end
end
