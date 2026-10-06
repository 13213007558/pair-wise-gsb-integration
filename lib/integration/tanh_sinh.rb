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

module Integration
  # Outcome of a tanh_sinh integration.
  #
  # [value]                 Best estimate of the integral.
  # [error]                 Worst-case error estimate:
  #                         max(discretization_error, tail_error).
  # [discretization_error]  Difference between the two finest step-size
  #                         levels (Float::INFINITY if only one level ran).
  # [tail_error]            Estimate of the truncated transformation tail:
  #                         magnitude of the first neglected term
  #                         (Float::INFINITY if the tail was never truncated).
  # [evaluations]           Number of integrand evaluations performed.
  # [levels]                Step-size halving levels completed.
  # [stop_reason]           Why the iteration stopped:
  #                         :converged, :max_levels, :max_evaluations or
  #                         :non_finite.
  TanhSinhResult = Struct.new(:value, :error, :discretization_error,
                              :tail_error, :evaluations, :levels,
                              :stop_reason) do
    # True only when both the discretization error and the transformation
    # tail truncation are within the requested tolerance.
    def converged?
      stop_reason == :converged
    end
  end

  class << self
    # Tanh-sinh (double exponential) quadrature over the finite interval
    # [t1, t2]. Handles integrable endpoint singularities such as log(x) or
    # 1/sqrt(x): the integrand is never evaluated at the endpoints, so the
    # caller does not need to provide finite function values there.
    #
    # Returns a TanhSinhResult reporting the estimate, the number of
    # integrand evaluations and the stop reason. A genuinely divergent
    # integral is never reported as converged.
    #
    # Options are:
    # [:tolerance]        Target absolute error for both the discretization
    #                     estimate and the tail truncation. Default: 1e-10.
    # [:max_levels]       Maximum number of step-size halvings. Default: 10.
    # [:max_evaluations]  Budget of integrand evaluations. Default: 100_000.
    def tanh_sinh(t1, t2, options = {}, &f)
      fail 'No function passed' unless block_given?
      fail 'Non-numeric bounds' unless t1.is_a?(Numeric) && t2.is_a?(Numeric)
      inf = infinite?(t1) || infinite?(t2) ||
            (t1.respond_to?(:infinite?) && t1.infinite?) ||
            (t2.respond_to?(:infinite?) && t2.infinite?)
      fail 'tanh_sinh does not support infinite bounds' if inf

      tolerance       = (options[:tolerance] || 1e-10).to_f
      max_levels      = options[:max_levels] || 10
      max_evaluations = options[:max_evaluations] || 100_000

      sign = 1
      lower = t1
      upper = t2
      if upper < lower
        lower, upper = upper, lower
        sign = -1
      end

      if lower == upper
        return TanhSinhResult.new(0.0, 0.0, 0.0, 0.0, 0, 0, :converged)
      end

      half = (upper - lower) / 2.0

      evaluations = 0
      tail_error  = Float::INFINITY

      # Transformed integrand f(x(t)) * x'(t) with
      # x(t) = center + half * tanh(pi/2 * sinh(t)).
      # The abscissa and the sech^2 weight factor are computed through
      # exp(-2 * |u|) so that x approaches the endpoints without
      # catastrophic cancellation. Points whose abscissa rounds onto an
      # endpoint carry a negligibly small weight and are skipped, so the
      # integrand is never called at the endpoints.
      # Aborts via throw(:abort, reason) on budget exhaustion or non-finite
      # values.
      sample = lambda do |t|
        u = (Math::PI / 2) * Math.sinh(t)
        eu = Math.exp(-2.0 * u.abs)
        ratio = eu / (1.0 + eu)
        x = if u <= 0
              lower + (upper - lower) * ratio
            else
              upper - (upper - lower) * ratio
            end
        next 0.0 if x == lower || x == upper
        sech2 = 4.0 * eu / ((1.0 + eu) * (1.0 + eu))
        weight = half * (Math::PI / 2) * Math.cosh(t) * sech2
        throw :abort, :max_evaluations if evaluations >= max_evaluations
        evaluations += 1
        term = f.call(x) * weight
        throw :abort, :non_finite unless term.is_a?(Numeric) && term.finite?
        term
      end

      # Sum sample.call(k * h) over both sides, truncating once a term
      # falls below the tolerance. With include_all, every integer k is
      # used (initial level); otherwise only odd k (refinement levels,
      # whose even multiples are reused from the previous level).
      # Updates tail_error with the magnitude of the first neglected term.
      sum_terms = lambda do |h, include_all|
        total = include_all ? sample.call(0.0) : 0.0
        increment = include_all ? 1 : 2
        level_tail = 0.0
        [1, -1].each do |side|
          k = 1
          loop do
            term = sample.call(side * k * h)
            if term.abs <= tolerance
              level_tail = [level_tail, term.abs].max
              break
            end
            total += term
            k += increment
          end
        end
        tail_error = level_tail
        total
      end

      h = 1.0
      sum = nil
      previous_sum = nil
      last_complete = nil
      level = 0
      discretization_error = Float::INFINITY

      aborted = catch(:abort) do
        loop do
          sum = h * sum_terms.call(h, level.zero?)
          sum += previous_sum / 2 unless previous_sum.nil?
          last_complete = sum

          if level >= 1
            discretization_error = (sum - previous_sum).abs
            if discretization_error <= tolerance && tail_error <= tolerance
              break :converged
            end
          end

          break :max_levels if level >= max_levels

          previous_sum = sum
          h /= 2
          level += 1
        end
      end

      stop_reason = aborted.is_a?(Symbol) ? aborted : :max_levels
      # An aborted level never truncated its tail; report it as unknown
      # rather than reuse the previous level's estimate.
      tail_error = Float::INFINITY if aborted == :max_evaluations || aborted == :non_finite
      value = last_complete || sum || 0.0
      error = [discretization_error, tail_error].max

      TanhSinhResult.new(sign * value, error, discretization_error,
                         tail_error, evaluations, level, stop_reason)
    end
  end
end
