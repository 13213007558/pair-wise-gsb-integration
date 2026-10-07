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
  # Result returned by the pure Ruby infinite interval integrator. The value
  # combines a finite numerical estimate with an independently computed interval
  # error and the caller supplied bound for the unintegrated tail.
  class InfiniteIntegralResult
    FIELDS = {
      value: :value,
      interval_error: :interval_error,
      computed_interval_error: :interval_error,
      tail_bound: :tail_bound,
      tail_upper_bound: :tail_bound,
      remainder_bound: :tail_bound,
      unproven_part: :unproven_part,
      unproven_tail: :unproven_part,
      interval: :interval,
      finite_interval: :interval,
      bounds: :interval,
      cutoffs: :cutoffs,
      cutoff: :cutoff,
      tolerance: :tolerance,
      tails: :tails
    }

    attr_reader :value, :interval_error, :tail_bound, :unproven_part,
                :interval, :cutoff, :cutoffs, :tolerance, :tails

    def initialize(attributes)
      @value = attributes.fetch(:value)
      @interval_error = attributes.fetch(:interval_error)
      @tail_bound = attributes[:tail_bound]
      @unproven_part = attributes.fetch(:unproven_part)
      @interval = attributes.fetch(:interval)
      @cutoff = attributes[:cutoff]
      @cutoffs = attributes[:cutoffs] || @cutoff
      @tolerance = attributes.fetch(:tolerance)
      @tails = attributes[:tails] || []
      @error_bound = finite_number?(@tail_bound) ? @interval_error + @tail_bound : Float::INFINITY
      @proven = finite_number?(@tail_bound) && @error_bound <= @tolerance && @unproven_part == 0.0
    end

    def proven?
      @proven
    end

    def unproven?
      !proven?
    end

    def error_bound
      @error_bound
    end

    def [](key)
      field = FIELDS[key] || (key.respond_to?(:to_sym) ? FIELDS[key.to_sym] : nil)
      public_send(field) if field
    rescue NoMethodError
      nil
    end

    def to_h
      {
        value: value,
        interval_error: interval_error,
        tail_bound: tail_bound,
        unproven_part: unproven_part,
        interval: interval,
        cutoff: cutoff,
        cutoffs: cutoffs,
        tolerance: tolerance,
        tails: tails,
        error_bound: error_bound,
        proven: proven?
      }
    end

    def to_f
      value.to_f
    end

    alias_method :computed_interval_error, :interval_error
    alias_method :tail_upper_bound, :tail_bound
    alias_method :remainder_bound, :tail_bound
    alias_method :unproven_tail, :unproven_part
    alias_method :unproven, :unproven_part
    alias_method :finite_interval, :interval

    private

    def finite_number?(value)
      value.is_a?(Numeric) && value.finite? && !value.nan?
    end
  end

  # Minus Infinity
  MInfinity = :minfinity

  # Infinity
  Infinity = :infinity

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
      value == Integration::Infinity || value == Integration::MInfinity ||
        (value.respond_to?(:infinite?) && value.infinite?)
    end

    # Get the integral for a function +f+, with bounds +t1+ and +t2+ given a
    # hash of +options+. If Ruby/GSL is available, you can use
    # +Integration::Minfinity+ and +Integration::Infinity+ as bounds. Method
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
    def integrate(t1, t2, options = {}, &f)
      fail 'No function passed' unless block_given?
      t1, t2, direction = normalize_bounds(t1, t2)
      inf_bounds = !direction.empty?

      fail 'Non-numeric bounds' unless ((t1.is_a? Numeric) && (t2.is_a? Numeric)) || inf_bounds

      unless inf_bounds
        lower_bound = [t1, t2].min
        upper_bound = [t1, t2].max
      end

      has_gsl = Integration.has_gsl?
      default_method = has_gsl ? :qag : :simpson
      options = { tolerance: 1e-10, initial_step: 16, step: 16, method: default_method }.merge(options)

      if inf_bounds && (!has_gsl || RUBY_METHODS.include?(options[:method]))
        unless options.key?(:tail_bound) || options.key?(:tail_upper_bound) ||
               options.key?(:remainder_bound) || options.key?(:upper_tail_bound) ||
               options.key?(:lower_tail_bound)
          fail 'A tail bound is required to prove an infinite Ruby integral'
        end

        result = integrate_infinite(t1, t2, options, &f)
        return result.value if result.proven?
        fail "Infinite integral was not proven: unproven_part=#{result.unproven_part}"
      elsif RUBY_METHODS.include? options[:method]
        integrate_ruby(lower_bound, upper_bound, options, &f)
      elsif GSL_METHODS.include? options[:method]
        if inf_bounds
          integrate_gsl(t1, t2, options, &f)
        else
          integrate_gsl(lower_bound, upper_bound, options, &f)
        end
      else
        fail "Unknown integration method \"#{options[:method]}\""
      end
    rescue RuntimeError => error
      raise unless error.message == 'Invalid infinite bounds'
      fail 'Non-numeric bounds'
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

    # Integrate an interval containing +MInfinity+ or +Infinity+ using pure Ruby.
    # When a tail bound callable is supplied it receives the finite cutoff and,
    # when it accepts two arguments, the direction +:upper+ or +:lower+. Returning
    # +nil+ or an infinite value marks the remaining tail as unproven rather than
    # treating a large finite endpoint as the complete integral. Without a tail
    # bound, only the finite computed interval and its error are reported.
    #
    # Options include +tolerance+, +method+, +initial_cutoff+, +cutoff+,
    # +extension_factor+, +max_extensions+, +split_at+, +tail_bound+,
    # +upper_tail_bound+, and +lower_tail_bound+.
    def integrate_infinite(t1, t2, options = {}, &f)
      fail 'No function passed' unless block_given?

      lower_bound, upper_bound, directions = normalize_bounds(t1, t2)
      fail 'Finite bounds are not supported by integrate_infinite' if directions.empty?
      fail 'Both bounds cannot be the same infinity' if directions.size == 2 && lower_bound == upper_bound

      options = {
        tolerance: 1e-10,
        initial_step: 16,
        step: 16,
        method: :adaptive_quadrature,
        initial_cutoff: nil,
        extension_factor: 2.0,
        max_extensions: 256,
        max_refinements: 100
      }.merge(options)

      fail 'Monte Carlo cannot provide a proven interval error' if options[:method] == :monte_carlo
      fail "Unknown integration method \"#{options[:method]}\"" unless RUBY_METHODS.include?(options[:method])

      split_point = options[:split_at]
      if directions.size == 2
        split_point = split_point.to_f if split_point
        fail 'split_at must be a finite number for doubly infinite integration' if split_point.nil? || !split_point.finite?
        ray_options = options.merge(error_budget: options[:tolerance] / 2.0)

        rays = [integrate_infinite_ray(MInfinity, split_point, ray_options, :lower, &f),
                integrate_infinite_ray(split_point, Infinity, ray_options, :upper, &f)]
        result = combine_infinite_rays([[MInfinity, split_point], [split_point, Infinity]], rays, options[:tolerance])
      else
        direction = directions.first
        ray_options = options.merge(error_budget: options[:tolerance])
        rays = [integrate_infinite_ray(lower_bound, upper_bound, ray_options, direction, &f)]
        result = combine_infinite_rays([[lower_bound, upper_bound]], rays, options[:tolerance])
      end

      result
    end

    alias_method :integrate_infinity, :integrate_infinite

    def integrate_ruby_with_error(lower_bound, upper_bound, options, &f)
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

      case method
      when :adaptive_quadrature
        max_level = options.key?(:max_level) ? options[:max_level] : 50
        Integration.adaptive_quadrature_with_error(lower_bound, upper_bound, tolerance, max_level, &f)
      when :romberg
        Integration.romberg_with_error(lower_bound, upper_bound, tolerance, &f)
      when :gauss
        n = initial_step > 10 ? 10 : initial_step
        value = method_obj.call(lower_bound, upper_bound, n, &f)
        return [value, 0.0] if n >= 10
        refined_n = n + 1
        refined = method_obj.call(lower_bound, upper_bound, refined_n, &f)
        [refined, (refined - value).abs]
      when :gauss_kronrod
        selected_points = points || 15
        value = method_obj.call(lower_bound, upper_bound, 10, selected_points == 15 ? nil : selected_points, &f)
        refined_points = selected_points >= 61 ? 61 : selected_points + 6
        refined = method_obj.call(lower_bound, upper_bound, 10, refined_points == 15 ? nil : refined_points, &f)
        [refined, (refined - value).abs]
      else
        current_step = initial_step
        value = method_obj.call(lower_bound, upper_bound, current_step, &f)
        maximum_step = options.key?(:max_refinements) ? initial_step + (step * options[:max_refinements]) : nil

        begin
          previous = value
          current_step += step
          value = method_obj.call(lower_bound, upper_bound, current_step, &f)
          error = (value - previous).abs
        end while error > tolerance && (maximum_step.nil? || current_step < maximum_step)

        [value, error]
      end
    end

    def integrate_ruby(lower_bound, upper_bound, options, &f)
      integrate_ruby_with_error(lower_bound, upper_bound, options, &f).first
    end

    private

    def normalize_bounds(t1, t2)
      lower = t1.is_a?(Float) && t1.infinite? ? (t1 < 0 ? MInfinity : Infinity) : t1
      upper = t2.is_a?(Float) && t2.infinite? ? (t2 < 0 ? MInfinity : Infinity) : t2
      directions = []
      directions << :lower if lower == MInfinity
      directions << :upper if upper == Infinity
      if lower == Infinity || upper == MInfinity ||
         (directions.include?(:lower) && directions.include?(:upper) && lower == upper) ||
         (lower == MInfinity && upper == MInfinity) || (lower == Infinity && upper == Infinity)
        fail 'Invalid infinite bounds'
      end
      fail 'Non-numeric bounds' unless valid_infinite_bounds?(lower, upper, directions)
      [lower, upper, directions]
    end

    def valid_infinite_bounds?(lower, upper, directions)
      return (lower.is_a?(Numeric) && upper.is_a?(Numeric)) if directions.empty?

      if directions == [:lower]
        upper.is_a?(Numeric) && upper.finite?
      elsif directions == [:upper]
        lower.is_a?(Numeric) && lower.finite?
      else
        lower == MInfinity && upper == Infinity
      end
    end

    def integrate_infinite_ray(lower_bound, upper_bound, options, direction, &f)
      anchor = direction == :upper ? lower_bound : upper_bound
      cutoff = initial_cutoff(options, anchor, direction)
      factor = options[:extension_factor].to_f.abs
      fail 'extension_factor must be greater than 1' unless factor > 1.0

      previous_cutoff = anchor
      total_value = 0.0
      total_error = 0.0
      interval = []
      segments = 0
      tail_bound = nil
      unproven_part = Float::INFINITY
      tail_function = tail_bound_function(options, direction)
      error_budget = options[:error_budget] || options[:tolerance]

      (options[:max_extensions] + 1).times do |extension|
        segment_start, segment_end = direction == :upper ? [previous_cutoff, cutoff] : [cutoff, previous_cutoff]
        segment_options = options.merge(tolerance: error_budget * ((factor - 1.0) / factor) / factor**segments)
        segment_value, segment_error = integrate_ruby_with_error(segment_start, segment_end, segment_options, &f)

        total_value += segment_value
        total_error += segment_error
        interval = ray_interval(lower_bound, upper_bound, cutoff, direction)
        tail_bound = evaluate_tail_bound(tail_function, cutoff, direction)
        segments += 1

        if tail_bound.nil?
          unproven_part = Float::INFINITY
          break
        elsif tail_bound <= error_budget + Float::EPSILON * error_budget - total_error
          unproven_part = 0.0
          break
        else
          unproven_part = tail_bound
        end

        if extension >= options[:max_extensions]
          break
        end

        previous_cutoff = cutoff
        cutoff = extend_cutoff(anchor, cutoff, factor, direction)
        if cutoff.nil? || !cutoff.finite?
          tail_bound = Float::INFINITY
          unproven_part = Float::INFINITY
          break
        end
      end

      tail_bound ||= Float::INFINITY
      if tail_bound.finite? && tail_bound <= error_budget + Float::EPSILON * error_budget - total_error
        unproven_part = 0.0
      elsif tail_bound.finite? && unproven_part != 0.0
        unproven_part = tail_bound
      end

      {
        direction: direction,
        value: total_value,
        interval_error: total_error,
        tail_bound: tail_bound,
        unproven_part: unproven_part,
        interval: interval,
        cutoff: cutoff,
        segments: segments
      }
    end

    def initial_cutoff(options, anchor, direction)
      configured = options[:initial_cutoff] || options[:cutoff]
      cutoff = configured ? configured.to_f : (anchor + (direction == :upper ? 1.0 : -1.0)).to_f
      fail 'initial_cutoff must be finite' unless cutoff.finite?
      unless direction == :upper ? cutoff > anchor : cutoff < anchor
        fail 'initial_cutoff must lie inside the infinite interval'
      end
      cutoff
    end

    def ray_interval(lower_bound, upper_bound, cutoff, direction)
      direction == :upper ? [lower_bound, cutoff] : [cutoff, upper_bound]
    end

    def extend_cutoff(anchor, cutoff, factor, direction)
      distance = (cutoff - anchor).abs * factor
      next_cutoff = direction == :upper ? anchor + distance : anchor - distance
      next_cutoff.finite? ? next_cutoff : nil
    end

    def tail_bound_function(options, direction)
      specific = direction == :upper ? options[:upper_tail_bound] : options[:lower_tail_bound]
      specific || options[:tail_bound] || options[:tail_upper_bound] || options[:remainder_bound]
    end

    def evaluate_tail_bound(function, cutoff, direction)
      return function.to_f if function.is_a?(Numeric) && function.finite?
      return Float::INFINITY if function.is_a?(Numeric) && function.infinite?
      return nil unless function.respond_to?(:call)

      arguments = function.arity == 2 ? [cutoff, direction] : [cutoff]
      bound = function.call(*arguments)
      return nil if bound.nil? || bound == false || !bound.is_a?(Numeric) || !bound.finite?
      fail 'Tail bounds must be non-negative' if bound < 0
      bound.to_f
    end

    def combine_infinite_rays(bounds, rays, tolerance)
      interval = [rays.first[:interval].first, rays.last[:interval].last]
      tail_bound = rays.inject(0.0) { |sum, ray| sum + ray[:tail_bound] }
      interval_error = rays.inject(0.0) { |sum, ray| sum + ray[:interval_error] }
      unproven_part = rays.all? { |ray| ray[:unproven_part] == 0.0 } ? 0.0 :
                      rays.any? { |ray| ray[:unproven_part] == Float::INFINITY } ? Float::INFINITY :
                      rays.inject(0.0) { |sum, ray| sum + ray[:unproven_part] }
      cutoffs = rays.size == 1 ? rays.first[:cutoff] : rays.map { |ray| ray[:cutoff] }

      InfiniteIntegralResult.new(
        value: rays.inject(0.0) { |sum, ray| sum + ray[:value] },
        interval_error: interval_error,
        tail_bound: tail_bound.finite? ? tail_bound : Float::INFINITY,
        unproven_part: unproven_part,
        interval: interval,
        cutoff: cutoffs,
        cutoffs: cutoffs,
        tolerance: tolerance,
        tails: rays
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
    public :has_gsl?
  end
end
