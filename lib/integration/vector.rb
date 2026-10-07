module Integration
  class ConvergenceError < StandardError
    attr_accessor :result
  end

  VectorComponentResult = Struct.new(:index, :value, :error, :abs_tolerance,
                                    :relative_tolerance, :target_error,
                                    :converged) do
    def converged?
      converged
    end
  end

  VectorIntegrationResult = Struct.new(:values, :errors, :components,
                                       :callback_count, :subdivisions,
                                       :intervals, :max_depth, :converged) do
    def converged?
      converged
    end

    def to_a
      values
    end
  end

  VectorInterval = Struct.new(:a, :b, :m, :fa, :fm, :fb, :value, :error,
                              :depth, :splittable)

  class << self
    # Adaptively integrate a callback that returns all components at once.
    #
    # Abs_tolerance and relative_tolerance may each be a scalar or an array.
    # Array tolerances are assigned by component position. Every component uses
    # the same sampled abscissas and interval partition, while convergence is
    # checked independently with abs_tolerance + relative_tolerance * |value|.
    def integrate_vector(t1, t2, options = {}, &f)
      fail ArgumentError, 'No function passed' unless block_given?
      fail ArgumentError, 'Non-numeric bounds' unless t1.is_a?(Numeric) && t2.is_a?(Numeric)
      fail ArgumentError, 'Vector integration does not support infinite bounds' if infinite?(t1) || infinite?(t2)
      unless t1.finite? && t2.finite?
        fail ArgumentError, 'Non-finite bounds'
      end

      lower_bound = [t1, t2].min.to_f
      upper_bound = [t1, t2].max.to_f
      orientation = t1 <= t2 ? 1.0 : -1.0
      abs_setting = tolerance_setting(options, :abs_tolerance, :tolerance, 1e-10)
      rel_setting = tolerance_setting(options, :relative_tolerance, :rel_tolerance, 0.0)
      max_subdivisions = Integer(options.fetch(:max_subdivisions, 1000))
      fail ArgumentError, 'max_subdivisions must be positive' unless max_subdivisions > 0

      dimension = nil
      abs_tolerances = nil
      relative_tolerances = nil
      values_at = {}
      callback_count = 0

      evaluate = lambda do |x|
        next values_at[x] if values_at.key?(x)

        callback_count += 1
        raw_values = f.call(x)
        unless raw_values.is_a?(Array)
          fail ArgumentError, "vector callback must return an Array, got #{raw_values.class} at x=#{x}"
        end

        values = raw_values.map do |raw_value|
          value = Float(raw_value)
          unless value.finite?
            fail ArgumentError, "vector callback returned non-finite value #{value.inspect} at x=#{x}"
          end
          value
        end

        if dimension.nil?
          dimension = values.length
          abs_tolerances = resolve_tolerance_setting(abs_setting, dimension, :abs_tolerance)
          relative_tolerances = resolve_tolerance_setting(rel_setting, dimension, :relative_tolerance)
        elsif values.length != dimension
          fail ArgumentError,
               "vector callback changed length from #{dimension} to #{values.length} at x=#{x}"
        end

        values_at[x] = values
      end

      first_values = evaluate.call(lower_bound)

      if dimension.zero?
        return vector_result([], [], [], callback_count, 0, 0, 0, true)
      end

      if lower_bound == upper_bound
        zero_values = Array.new(dimension, 0.0)
        zero_errors = Array.new(dimension, 0.0)
        components = vector_components(zero_values, zero_errors,
                                       abs_tolerances, relative_tolerances)
        return vector_result(zero_values, zero_errors, components,
                             callback_count, 0, 0, 0, true)
      end

      middle = (lower_bound + upper_bound) / 2.0
      left_middle = lower_bound + (middle - lower_bound) / 2.0
      right_middle = middle + (upper_bound - middle) / 2.0

      fa = first_values
      fm = evaluate.call(middle)
      fb = evaluate.call(upper_bound)
      f_left_middle = evaluate.call(left_middle)
      f_right_middle = evaluate.call(right_middle)

      width = upper_bound - lower_bound
      coarse = scale_vector(add_vectors(fa, scale_vector(fm, 4.0), fb), width / 6.0)
      left_value = scale_vector(add_vectors(fa, scale_vector(f_left_middle, 4.0), fm), width / 12.0)
      right_value = scale_vector(add_vectors(fm, scale_vector(f_right_middle, 4.0), fb), width / 12.0)
      fine = add_vectors(left_value, right_value)
      root_error = scale_vector(subtract_vectors(fine, coarse), 1.0 / 15.0).map(&:abs)

      values = fine
      errors = root_error
      intervals = [
        VectorInterval.new(lower_bound, middle, left_middle, fa, f_left_middle, fm,
                           left_value, scale_vector(root_error, 0.5), 1, true),
        VectorInterval.new(middle, upper_bound, right_middle, fm, f_right_middle, fb,
                           right_value, scale_vector(root_error, 0.5), 1, true)
      ]
      subdivisions = 1
      max_depth = 1

      loop do
        targets = component_targets(values, abs_tolerances, relative_tolerances)
        converged = errors.each_with_index.all? { |error, index| error <= targets[index] }
        result = nil

        if converged
          result_values = scale_vector(values, orientation)
          components = vector_components(values, errors, abs_tolerances,
                                         relative_tolerances, targets, result_values)
          return vector_result(result_values, errors, components, callback_count,
                               subdivisions, intervals.length, max_depth, true)
        end

        if subdivisions >= max_subdivisions
          result_values = scale_vector(values, orientation)
          components = vector_components(values, errors, abs_tolerances,
                                         relative_tolerances, targets, result_values)
          result = vector_result(result_values, errors, components, callback_count,
                                 subdivisions, intervals.length, max_depth, false)
          error = ConvergenceError.new('Vector integral did not converge')
          error.result = result
          fail error
        end

        interval = select_interval(intervals, errors, targets)
        unless interval
          result_values = scale_vector(values, orientation)
          components = vector_components(values, errors, abs_tolerances,
                                         relative_tolerances, targets, result_values)
          result = vector_result(result_values, errors, components, callback_count,
                                 subdivisions, intervals.length, max_depth, false)
          error = ConvergenceError.new('Vector integral reached machine precision')
          error.result = result
          fail error
        end

        split = split_interval(interval, evaluate)

        if split.fetch(:degenerate)
          errors = subtract_vectors(errors, interval.error)
          interval.error = split.fetch(:child_error)
          interval.splittable = false
          errors = add_vectors(errors, interval.error)
        else
          subdivisions += 1
          max_depth = [max_depth, interval.depth + 1].max
          replace_interval!(intervals, interval, split.fetch(:left), split.fetch(:right))

          refined_values = add_vectors(split.fetch(:left).value,
                                       split.fetch(:right).value)
          values = add_vectors(values, subtract_vectors(refined_values, interval.value))
          errors = subtract_vectors(errors, interval.error)
          errors = add_vectors(errors, split.fetch(:left).error)
          errors = add_vectors(errors, split.fetch(:right).error)
        end
      end
    end

    private

    def tolerance_setting(options, key, legacy_key, default)
      value = options.fetch(key, options.fetch(legacy_key, default))
      if value.is_a?(Array)
        value.map { |item| validate_tolerance(item, key) }
      else
        validate_tolerance(value, key)
      end
    end

    def validate_tolerance(value, key)
      tolerance = Float(value)
      unless tolerance.finite? && tolerance >= 0
        fail ArgumentError, "#{key} must be finite and non-negative"
      end
      tolerance
    end

    def resolve_tolerance_setting(setting, dimension, key)
      values = setting.is_a?(Array) ? setting : Array.new(dimension, setting)
      if values.length != dimension
        fail ArgumentError,
             "#{key} has #{values.length} entries but callback returned #{dimension} components"
      end
      values
    end

    def component_targets(values, abs_tolerances, relative_tolerances)
      values.each_with_index.map do |value, index|
        abs_tolerances[index] + relative_tolerances[index] * value.abs
      end
    end

    def select_interval(intervals, errors, targets)
      selected = nil
      best_score = -1.0

      intervals.each do |candidate|
        next unless candidate.splittable

        score = 0.0
        errors.each_with_index do |total_error, index|
          next unless total_error > targets[index]
          ratio = if targets[index] > 0
                    candidate.error[index] / targets[index]
                  else
                    candidate.error[index] > 0 ? Float::INFINITY : 0.0
                  end
          score = ratio if ratio > score
        end
        selected = candidate if score > best_score
        best_score = score if score > best_score
      end

      selected
    end

    def split_interval(interval, evaluate)
      left_midpoint = interval.a + (interval.m - interval.a) / 2.0
      right_midpoint = interval.m + (interval.b - interval.m) / 2.0
      degenerate = left_midpoint == interval.a || left_midpoint == interval.m ||
                   right_midpoint == interval.m || right_midpoint == interval.b
      f_left_midpoint = evaluate.call(left_midpoint)
      f_right_midpoint = evaluate.call(right_midpoint)
      width = interval.b - interval.a

      left_value = scale_vector(
        add_vectors(interval.fa, scale_vector(f_left_midpoint, 4.0), interval.fm),
        width / 12.0
      )
      right_value = scale_vector(
        add_vectors(interval.fm, scale_vector(f_right_midpoint, 4.0), interval.fb),
        width / 12.0
      )
      refined_value = add_vectors(left_value, right_value)
      refined_error = scale_vector(subtract_vectors(refined_value, interval.value), 1.0 / 15.0).map(&:abs)
      child_error = scale_vector(refined_error, 0.5)

      if degenerate
        return { degenerate: true, child_error: refined_error }
      end

      left_interval = VectorInterval.new(interval.a, interval.m, left_midpoint,
                                         interval.fa, f_left_midpoint, interval.fm,
                                         left_value, child_error, interval.depth + 1, true)
      right_interval = VectorInterval.new(interval.m, interval.b, right_midpoint,
                                          interval.fm, f_right_midpoint, interval.fb,
                                          right_value, child_error, interval.depth + 1, true)

      { degenerate: false, left: left_interval, right: right_interval }
    end

    def replace_interval!(intervals, interval, left_interval, right_interval)
      index = intervals.index(interval)
      intervals[index, 1] = [left_interval, right_interval]
    end

    def add_vectors(*vectors)
      vectors.each_with_object([]) do |vector, result|
        vector.each_with_index { |value, index| result[index] = result[index].to_f + value }
      end
    end

    def subtract_vectors(left, right)
      left.each_with_index.map { |value, index| value - right[index] }
    end

    def scale_vector(vector, factor)
      vector.map { |value| value * factor }
    end

    def vector_components(values, errors, abs_tolerances, relative_tolerances,
                          targets = component_targets(values, abs_tolerances, relative_tolerances),
                          result_values = values)
      result_values.each_with_index.map do |value, index|
        VectorComponentResult.new(index, value, errors[index], abs_tolerances[index],
                                  relative_tolerances[index], targets[index],
                                  errors[index] <= targets[index])
      end
    end

    def vector_result(values, errors, components, callback_count, subdivisions,
                      interval_count, max_depth, converged)
      VectorIntegrationResult.new(values, errors, components, callback_count,
                                  subdivisions, interval_count, max_depth, converged)
    end
  end
end
