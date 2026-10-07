# Copyright (c) 2005  Beng (original code)
#               2011  Claudio Bustos
#               XXXX -> Add new developers
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE X CONSORTIUM BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
# OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Except as contained in this notice, the name of the Beng shall not
# be used in advertising or otherwise to promote the sale, use or other dealings
# in this Software without prior written authorization from Beng.

require 'integration/methods'

# Diverse integration methods
# Use Integration.integrate as wrapper to direct access to methods
#
module Integration

  # Minus Infinity
  MInfinity = :minfinity

  # Infinity
  Infinity = :infinity

  # Result of a pure-Ruby integration over an unbounded interval. The object
  # behaves as the computed finite-part value for simple numeric operations.
  class InfiniteIntegralResult < Numeric
    attr_reader :value, :interval_error, :tail_bound, :unproven,
                :lower_cutoff, :upper_cutoff, :iterations

    def initialize(value, interval_error, tail_bound, unproven, options)
      @value = value
      @interval_error = interval_error
      @tail_bound = tail_bound
      @unproven = unproven
      @lower_cutoff = options[:lower_cutoff]
      @upper_cutoff = options[:upper_cutoff]
      @iterations = options[:iterations]
    end

    def proven?
      unproven == 0
    end

    def to_f
      value.to_f
    end

    def to_i
      to_f.to_i
    end

    def coerce(other)
      [other, to_f]
    end

    def <=>(other)
      to_f <=> other
    end

    def ==(other)
      other == to_f
    end

    def +(other)
      to_f + other
    end

    def -(other)
      to_f - other
    end

    def *(other)
      to_f * other
    end

    def /(other)
      to_f / other
    end

    def abs
      to_f.abs
    end

    def to_h
      {
        value: value,
        interval_error: interval_error,
        tail_bound: tail_bound,
        unproven: unproven,
        lower_cutoff: lower_cutoff,
        upper_cutoff: upper_cutoff,
        iterations: iterations
      }
    end

    def inspect
      "#<Integration::InfiniteIntegralResult value=#{value} " \
      "interval_error=#{interval_error} tail_bound=#{tail_bound} " \
      "unproven=#{unproven}>"
    end
    alias_method :to_s, :inspect
  end

  # Pure Ruby methods available.
  RUBY_METHODS = [:rectangle, :trapezoid, :simpson, :adaptive_quadrature,
                 :gauss, :romberg, :monte_carlo, :gauss_kronrod,
                 :simpson3by8, :boole, :open_trapezoid, :milne]

  # Methods available when using the `rb-gsl` gem.
  GSL_METHODS = [:qng, :qag]

  class << self

    # Check if `value` is plus or minus infinity.
    #
    # @param value Value to be tested.
    def infinite?(value)
      value == Integration::Infinity || value == Integration::MInfinity
    end

    # Get the integral for a function +f+, with bounds +t1+ and +t2+ given a
    # hash of +options+. Infinite bounds use GSL when available. Without GSL,
    # they use {#integrate_infinite_ruby}, which can report a tail bound.
    #
    # Options are:
    # [:tolerance]    Maximum difference between real and calculated integral.
    #                 Default: 1e-10.
    # [:initial_step] Initial number of subdivisions.
    # [:step]         Subdivition increment on each iteration.
    # [:method]       Integration method.
    #
    # Available methods are:
    #
    # [:rectangle] for [:initial_step+:step*iteration] quadrilateral subdivisions.
    # [:trapezoid] for [:initial_step+:step*iteration] trapezoid-al subdivisions.
    # [:simpson]   for [:initial_step+:step*iteration] parabolic subdivisions.
    # [:adaptive_quadrature] for recursive appoximations until error [tolerance].
    # [:gauss] [:initial_step+:step*iteration] weighted subdivisons using
    # translated -1 -> +1 endpoints.
    # [:romberg] extrapolation of recursion approximation until error < [tolerance].
    # [:monte_carlo] make [:initial_step+:step*iteration] random samples, and
    # check for above/below curve.
    # [:qng] GSL QNG non-adaptive Gauss-Kronrod integration.
    # [:qag] GSL QAG adaptive integration, with support for infinite bounds.
    #
    # Infinite-bound options for the pure-Ruby implementation are documented on
    # {#integrate_infinite_ruby}.
    def integrate(t1, t2, options = {}, &f)
      inf_bounds = (infinite?(t1) || infinite?(t2))

      fail 'No function passed' unless block_given?
      fail 'Non-numeric bounds' unless ((t1.is_a? Numeric) && (t2.is_a? Numeric)) || inf_bounds

      if inf_bounds
        lower_bound = t1
        upper_bound = t2
        options[:method] = :qag if options[:method].nil? && has_gsl?
      else
        lower_bound = [t1, t2].min
        upper_bound = [t1, t2].max
      end

      def_method = (Integration.has_gsl?) ? :qag : :simpson
      default_opts = { tolerance: 1e-10, initial_step: 16, step: 16, method: def_method }
      options = default_opts.merge(options)

      if inf_bounds && use_ruby_infinite_integration?(options)
        integrate_infinite_ruby(lower_bound, upper_bound, options, &f)
      elsif RUBY_METHODS.include? options[:method]
        integrate_ruby(lower_bound, upper_bound, options, &f)
      elsif GSL_METHODS.include? options[:method]
        integrate_gsl(lower_bound, upper_bound, options, &f)
      else
        fail "Unknown integration method \"#{options[:method]}\""
      end
    end

    # Pure-Ruby integration on an unbounded interval.
    #
    # The interval is expanded from +cutoff+ until the sum of the computed
    # finite interval error and the caller-supplied tail bound is no larger
    # than +tolerance+. A tail bound is a callable that receives a cutoff point
    # and returns a nonnegative upper bound on the absolute value of the
    # remaining integral. Use +lower_tail_bound+ and +upper_tail_bound+ for a
    # doubly infinite interval. A single +tail_bound+ may accept an optional
    # second +:lower+ or +:upper+ direction argument.
    #
    # If no finite tail bound is available, the returned result is never marked
    # as proven; its +unproven+ value is Infinity.
    def integrate_infinite_ruby(lower_bound, upper_bound, options = {}, &f)
      fail 'No function passed' unless block_given?

      tolerance = options.fetch(:tolerance, 1e-10).to_f
      factor = options.fetch(:cutoff_factor, 10.0).to_f
      max_iterations = options.fetch(:max_iterations, 32).to_i
      max_depth = options.fetch(:max_depth, 100).to_i
      fail 'tolerance must be positive' unless tolerance > 0
      fail 'cutoff_factor must be greater than one' unless factor > 1
      fail 'max_iterations must be positive' unless max_iterations > 0

      lower_infinite = infinite?(lower_bound)
      upper_infinite = infinite?(upper_bound)
      fail 'At least one bound must be infinite' unless lower_infinite || upper_infinite
      fail 'QNG does not support infinite bounds' if options[:method] == :qng
      unless lower_infinite || lower_bound.is_a?(Numeric)
        fail 'Non-numeric lower bound'
      end
      unless upper_infinite || upper_bound.is_a?(Numeric)
        fail 'Non-numeric upper bound'
      end
      if !lower_infinite && !upper_infinite && lower_bound > upper_bound
        lower_bound, upper_bound = upper_bound, lower_bound
      end

      if lower_infinite && upper_infinite
        center = options.fetch(:center, 0.0).to_f
        initial = options.fetch(:cutoff, [center.abs + 1.0, 1.0].max).to_f
        current_lower = center - initial
        current_upper = center + initial
      elsif upper_infinite
        initial = options.fetch(:cutoff, lower_bound.to_f + 1.0).to_f
        fail 'cutoff must be greater than lower bound' unless initial > lower_bound
        current_lower = lower_bound.to_f
        current_upper = initial
        initial = current_upper - current_lower
      else
        initial = options.fetch(:cutoff, upper_bound.to_f - 1.0).to_f
        fail 'cutoff must be less than upper bound' unless initial < upper_bound
        current_lower = initial
        current_upper = upper_bound.to_f
        initial = current_upper - current_lower
      end

      lower_bound_fn = options[:lower_tail_bound]
      upper_bound_fn = options[:upper_tail_bound]
      common_bound_fn = options[:tail_bound]
      if lower_infinite && lower_bound_fn.nil? && common_bound_fn
        lower_bound_fn = wrap_common_tail_bound(common_bound_fn, :lower)
      end
      if upper_infinite && upper_bound_fn.nil? && common_bound_fn
        upper_bound_fn = wrap_common_tail_bound(common_bound_fn, :upper)
      end

      value = 0.0
      interval_error = 0.0
      segment_tolerance = tolerance / (2.0 * max_iterations)
      previous_distance = initial
      tail_error = Float::INFINITY
      iterations = 0
      distance = initial

      max_iterations.times do |iteration|
        iterations = iteration + 1

        if iteration == 0
          if lower_infinite && upper_infinite
            add_segment(current_lower, center, segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value += segment_value
              interval_error += segment_error
            end
            add_segment(center, current_upper, segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value += segment_value
              interval_error += segment_error
            end
          elsif upper_infinite
            add_segment(current_lower, current_upper, segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value = segment_value
              interval_error = segment_error
            end
          else
            add_segment(current_lower, current_upper, segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value = segment_value
              interval_error = segment_error
            end
          end
        else
          distance = initial * (factor ** iterations)
          u_min = initial / distance
          u_max = iteration == 1 ? 1.0 : initial / previous_distance
          current_lower = center - distance if lower_infinite && upper_infinite
          current_upper = center + distance if lower_infinite && upper_infinite
          current_upper = lower_bound.to_f + distance if upper_infinite && !lower_infinite
          current_lower = upper_bound.to_f - distance if lower_infinite && !upper_infinite

          if lower_infinite
            finite_point = upper_infinite ? center : upper_bound.to_f
            add_distance_segment(finite_point, -1, initial, u_min, u_max,
                                 segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value += segment_value
              interval_error += segment_error
            end
          end
          if upper_infinite
            finite_point = lower_infinite ? center : lower_bound.to_f
            add_distance_segment(finite_point, 1, initial, u_min, u_max,
                                 segment_tolerance, max_depth, f) do |segment_value, segment_error|
              value += segment_value
              interval_error += segment_error
            end
          end
        end

        distance = initial if iteration == 0

        tail_error = tail_value(lower_infinite, current_lower, lower_bound_fn,
                                upper_infinite, current_upper, upper_bound_fn)

        if finite_error?(interval_error) && finite_error?(tail_error) &&
           interval_error + tail_error <= tolerance
          return infinite_result(value, interval_error, tail_error, 0,
                                 current_lower, current_upper, iterations)
        end

        break unless finite_error?(interval_error) && finite_error?(tail_error)
        break if iterations == max_iterations

        previous_distance = distance
      end

      unproven = if finite_error?(interval_error) && finite_error?(tail_error)
        [interval_error + tail_error - tolerance, 0.0].max
      else
        Float::INFINITY
      end
      infinite_result(value, interval_error, tail_error, unproven,
                      current_lower, current_upper, iterations)
    end

    # Public entry point for the pure-Ruby infinite-bound implementation.
    def integrate_infinite(lower_bound, upper_bound, options = {}, &f)
      integrate_infinite_ruby(lower_bound, upper_bound, options, &f)
    end

    # Integrate using the GSL bindings.
    def integrate_gsl(lower_bound, upper_bound, options, &f)
      f = GSL::Function.alloc(&f)
      method = options[:method]
      tolerance = options[:tolerance]

      if (method == :qag)
        w = GSL::Integration::Workspace.alloc

        val = if infinite?(lower_bound) && infinite?(upper_bound)
          f.qagi([tolerance, 0.0], 1000, w)
        elsif infinite?(lower_bound)
          f.qagil(upper_bound, [tolerance, 0], w)
        elsif infinite?(upper_bound)
          f.qagiu(lower_bound, [tolerance, 0], w)
        else
          f.qag([lower_bound, upper_bound], [tolerance, 0.0], GSL::Integration::GAUSS61, w)
        end

      elsif (method == :qng)
        val = f.qng([lower_bound, upper_bound], [tolerance, 0.0])

      else
        fail "Unknown integration method \"#{method}\""
      end

      val[0]
    end

    def integrate_ruby(lower_bound, upper_bound, options, &f)
      method = options[:method]
      tolerance = options[:tolerance]
      initial_step = options[:initial_step]
      step = options[:step]
      points = options[:points]

      begin
        method_obj = Integration.method(method.to_s.downcase)
      rescue
        raise "Unknown integration method \"#{method}\""
      end

      current_step = initial_step

      if [:adaptive_quadrature, :romberg, :gauss, :gauss_kronrod].include? method
        if (method == :gauss)
          initial_step = 10 if initial_step > 10
          tolerance = initial_step
          method_obj.call(lower_bound, upper_bound, tolerance, &f)
        elsif (method == :gauss_kronrod)
          initial_step = 10 if initial_step > 10
          tolerance = initial_step
          points = points unless points.nil?
          method_obj.call(lower_bound, upper_bound, tolerance, points, &f)
        else
          method_obj.call(lower_bound, upper_bound, tolerance, &f)
        end
      else
        value = method_obj.call(lower_bound, upper_bound, current_step, &f)
        previous = value + (tolerance * 2)
        diffs = []

        while (previous - value).abs > tolerance
          diffs.push((previous - value).abs)
          current_step += step
          previous = value
          value = method_obj.call(lower_bound, upper_bound, current_step, &f)
        end

        value
      end
    end

    def use_ruby_infinite_integration?(options)
      RUBY_METHODS.include?(options[:method]) ||
        options.key?(:tail_bound) ||
        options.key?(:lower_tail_bound) ||
        options.key?(:upper_tail_bound) ||
        !has_gsl?
    end

    def wrap_common_tail_bound(function, direction)
      if function.respond_to?(:arity) && function.arity.abs == 2
        proc { |cutoff| function.call(cutoff, direction) }
      else
        function
      end
    end

    def add_segment(lower_bound, upper_bound, tolerance, max_depth, function)
      segment_value, segment_error = adaptive_simpson_with_error(
        lower_bound, upper_bound, tolerance, max_depth, &function
      )
      yield(segment_value, segment_error)
    end

    def add_distance_segment(finite_point, direction, initial, s_min, s_max,
      tolerance, max_depth, function)
      transformed = proc do |s|
        x = finite_point + direction * (initial / s)
        function.call(x) * initial / (s * s)
      end
      segment_value, segment_error = adaptive_simpson_with_error(
        s_min, s_max, tolerance, max_depth, &transformed
      )
      yield(segment_value, segment_error)
    end

    def tail_value(lower_infinite, lower_cutoff, lower_function,
                   upper_infinite, upper_cutoff, upper_function)
      values = []
      if lower_infinite
        return Float::INFINITY unless lower_function
        values << validate_tail_bound(lower_function.call(lower_cutoff), lower_cutoff)
      end
      if upper_infinite
        return Float::INFINITY unless upper_function
        values << validate_tail_bound(upper_function.call(upper_cutoff), upper_cutoff)
      end
      values.inject(0.0) { |sum, value| sum + value }
    end

    def validate_tail_bound(value, cutoff)
      valid = value.is_a?(Numeric) &&
              (!value.respond_to?(:nan?) || !value.nan?) &&
              value >= 0
      fail "Tail bound at #{cutoff} must be a nonnegative number" unless valid
      value.to_f
    end

    def finite_error?(value)
      value.is_a?(Numeric) && value.finite? rescue false
    end

    def infinite_result(value, interval_error, tail_bound, unproven,
                        lower_cutoff, upper_cutoff, iterations)
      InfiniteIntegralResult.new(
        value,
        interval_error,
        tail_bound,
        unproven,
        lower_cutoff: lower_cutoff,
        upper_cutoff: upper_cutoff,
        iterations: iterations
      )
    end

    # Check if GSL is available. Require the library if it is present. Return a
    # boolean indicating its availability.
    #
    # @return [Boolean] Whether GSL is available.
    def has_gsl?
      gsl_available = '@@gsl'
      if class_variable_defined? gsl_available
        class_variable_get(gsl_available)
      else
        begin
          require 'gsl'
          class_variable_set(gsl_available, true)
        rescue LoadError
          class_variable_set(gsl_available, false)
        end
      end
    end
  end
end
