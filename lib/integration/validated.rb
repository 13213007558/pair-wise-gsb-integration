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

module Integration
  # Result of Integration.validated_integrate.
  #
  # Unlike the plain quadrature methods, a validated result is a *proof*:
  # the true value of the integral is guaranteed to lie between #lower and
  # #upper, provided the caller-supplied enclosure was valid for every
  # subinterval it was queried on.
  class ValidatedIntegral
    # Guaranteed lower bound of the integral.
    attr_reader :lower

    # Guaranteed upper bound of the integral.
    attr_reader :upper

    # Number of enclosure evaluations used.
    attr_reader :evaluations

    # :converged when the requested tolerance was met,
    # :budget_exhausted when max_evaluations ran out first (the returned
    # enclosure is still valid, just wider than requested).
    attr_reader :status

    def initialize(lower, upper, evaluations, status)
      @lower = lower
      @upper = upper
      @evaluations = evaluations
      @status = status
    end

    # Width of the enclosing interval (upper - lower).
    def width
      upper - lower
    end

    # True when the requested tolerance was reached.
    def converged?
      status == :converged
    end

    # True when +value+ lies inside the guaranteed enclosure.
    def include?(value)
      lower <= value && value <= upper
    end

    # Midpoint of the enclosure (best point estimate).
    def midpoint
      (lower + upper) / 2
    end

    def to_a
      [lower, upper]
    end
  end

  # Max-heap of subinterval contributions, keyed by the width each
  # subinterval contributes to the total enclosure. Used to greedily
  # refine the subinterval that currently dominates the error.
  # @private
  class ContributionHeap
    def initialize
      @entries = []
    end

    def empty?
      @entries.empty?
    end

    def push(entry)
      @entries << entry
      i = @entries.size - 1
      while i > 0
        parent = (i - 1) / 2
        break if @entries[parent][0] >= entry[0]

        @entries[i] = @entries[parent]
        i = parent
      end
      @entries[i] = entry
    end

    def pop
      top = @entries[0]
      last = @entries.pop
      unless @entries.empty?
        i = 0
        loop do
          left = (2 * i) + 1
          right = left + 1
          break if left >= @entries.size

          child = right < @entries.size && @entries[right][0] > @entries[left][0] ? right : left
          break if @entries[child][0] <= last[0]

          @entries[i] = @entries[child]
          i = child
        end
        @entries[i] = last
      end
      top
    end
  end

  class << self
    # Validated (guaranteed enclosure) integration.
    #
    # Computes rigorous lower and upper bounds for the integral of a
    # function over [t1, t2]. Instead of a pointwise function, the caller
    # supplies an *enclosure oracle*: a block receiving a subinterval
    # (+a+, +b+) and returning a pair [lo, hi] such that
    #
    #   lo <= f(x) <= hi   for all x in [a, b]
    #
    # The algorithm adaptively bisects the subinterval with the widest
    # contribution and accumulates the per-subinterval bounds, so the
    # returned ValidatedIntegral always contains the true integral.
    #
    # Sampling a plain Float callback (as the other methods do) cannot
    # prove anything about the sampled error; this method never pretends
    # otherwise. Bounds are either computed exactly (Rational/Integer
    # arithmetic) or with explicit outward rounding of every Float
    # operation, tracked through a running error bound.
    #
    # Options are:
    # [:tolerance]       Target for (upper - lower). Default: 1e-10.
    # [:max_evaluations] Budget of enclosure oracle calls. When the budget
    #                    runs out before the tolerance is met, the current
    #                    (still valid) enclosure is returned with status
    #                    :budget_exhausted. Default: 10_000.
    # [:arithmetic]      :auto (default), :exact or :float.
    #                    :exact requires Rational/Integer bounds and
    #                    enclosures and performs exact arithmetic; it
    #                    raises if the oracle ever returns a Float.
    #                    :float uses Float arithmetic with explicit
    #                    outward rounding. :auto picks :exact when the
    #                    bounds and the first enclosure are exact,
    #                    :float otherwise.
    #
    # Example: integral of x^2 over [0, 1] is 1/3, enclosed exactly.
    #
    #   result = Integration.validated_integrate(Rational(0), Rational(1),
    #                                            tolerance: Rational(1, 1000)) do |a, b|
    #     [a * a, b * b] # x^2 is increasing on [0, 1]
    #   end
    #   result.include?(Rational(1, 3)) # => true
    def validated_integrate(t1, t2, options = {}, &enclosure)
      fail 'No enclosure passed' unless block_given?
      fail 'Validated integration does not support infinite bounds' if infinite?(t1) || infinite?(t2)
      fail 'Non-numeric bounds' unless t1.is_a?(Numeric) && t2.is_a?(Numeric)

      lower_bound = [t1, t2].min
      upper_bound = [t1, t2].max

      tolerance = options.fetch(:tolerance, 1e-10)
      max_evaluations = options.fetch(:max_evaluations, 10_000)
      arithmetic = options.fetch(:arithmetic, :auto)

      fail 'tolerance must be non-negative' unless tolerance.is_a?(Numeric) && tolerance >= 0
      fail 'max_evaluations must be a positive integer' unless max_evaluations.is_a?(Integer) && max_evaluations >= 1
      fail "Unknown arithmetic mode \"#{arithmetic}\"" unless [:auto, :exact, :float].include?(arithmetic)

      first_lo, first_hi = validated_enclosure(enclosure, lower_bound, upper_bound)

      exact = case arithmetic
              when :exact then true
              when :float then false
              else [lower_bound, upper_bound, first_lo, first_hi].all? { |v| exact_value?(v) }
              end

      if exact
        # Bisection midpoints must stay exact; Integer division would not.
        lower_bound = lower_bound.to_r
        upper_bound = upper_bound.to_r
      else
        lower_bound = lower_bound.to_f
        upper_bound = upper_bound.to_f
      end

      state = exact ? ExactEnclosureState.new : FloatEnclosureState.new

      evaluations = 1
      heap = ContributionHeap.new
      push_contribution(heap, state, lower_bound, upper_bound, first_lo, first_hi, exact)

      status = state.width <= tolerance ? :converged : nil

      while status.nil?
        if evaluations + 2 > max_evaluations || heap.empty?
          status = :budget_exhausted
          break
        end

        _key, a, b, lo_term, lo_err, hi_term, hi_err = heap.pop
        mid = (a + b) / 2

        if mid <= a || mid >= b
          # Subinterval cannot be bisected any further (adjacent floats);
          # its contribution stays in the totals and refinement stops.
          status = :budget_exhausted
          break
        end

        state.remove(lo_term, lo_err, hi_term, hi_err)

        left_lo, left_hi = validated_enclosure(enclosure, a, mid)
        right_lo, right_hi = validated_enclosure(enclosure, mid, b)
        evaluations += 2

        push_contribution(heap, state, a, mid, left_lo, left_hi, exact)
        push_contribution(heap, state, mid, b, right_lo, right_hi, exact)

        status = :converged if state.width <= tolerance
      end

      ValidatedIntegral.new(state.lower, state.upper, evaluations, status)
    end

    private

    def exact_value?(value)
      value.is_a?(Rational) || value.is_a?(Integer)
    end

    # Calls the enclosure oracle and validates its contract.
    def validated_enclosure(enclosure, a, b)
      result = enclosure.call(a, b)
      unless result.is_a?(Array) && result.size == 2
        fail "Enclosure must return a [lower, upper] pair, got #{result.inspect}"
      end

      lo, hi = result
      fail 'Enclosure bounds must be numeric' unless lo.is_a?(Numeric) && hi.is_a?(Numeric)
      fail "Invalid enclosure [#{lo}, #{hi}] on [#{a}, #{b}]: lower > upper" unless lo <= hi

      [lo, hi]
    end

    # Computes the contribution of one subinterval to the totals and
    # pushes it onto the refinement heap.
    def push_contribution(heap, state, a, b, lo, hi, exact)
      if exact && !(exact_value?(lo) && exact_value?(hi))
        fail 'Enclosure returned inexact (Float) bounds in exact arithmetic mode; ' \
             'use arithmetic: :float or return Rational bounds'
      end

      lo_term, lo_err, hi_term, hi_err = state.contribution(lo, hi, a, b)
      state.add(lo_term, lo_err, hi_term, hi_err)
      key = (hi_term + hi_err) - (lo_term - lo_err)
      heap.push([key, a, b, lo_term, lo_err, hi_term, hi_err])
    end
  end

  # Exact (Rational/Integer) accumulation of enclosure contributions.
  # Every operation is exact, so the totals are rigorous by construction.
  # @private
  class ExactEnclosureState
    def initialize
      @total_lo = 0
      @total_hi = 0
    end

    def contribution(lo, hi, a, b)
      w = b - a
      [lo * w, 0, hi * w, 0]
    end

    def add(lo_term, _lo_err, hi_term, _hi_err)
      @total_lo += lo_term
      @total_hi += hi_term
    end

    def remove(lo_term, _lo_err, hi_term, _hi_err)
      @total_lo -= lo_term
      @total_hi -= hi_term
    end

    def lower
      @total_lo
    end

    def upper
      @total_hi
    end

    def width
      @total_hi - @total_lo
    end
  end

  # Float accumulation with explicit outward rounding. Every operation
  # enlarges a running error bound by the worst-case rounding error, so
  # [lower, upper] is guaranteed to contain the exact infinite-precision
  # sum of the enclosure contributions.
  # @private
  class FloatEnclosureState
    # Conservative unit roundoff used for all error bounds
    # (Float::EPSILON is twice the true unit roundoff for round-to-nearest).
    U = Float::EPSILON

    def initialize
      @sum_lo = 0.0
      @err_lo = 0.0
      @sum_hi = 0.0
      @err_hi = 0.0
    end

    # Contribution of one subinterval, with a rigorous error bound.
    # lo/hi may be any Numeric; non-Float values are converted with .to_f
    # and the conversion error is accounted for.
    def contribution(lo, hi, a, b)
      w = b - a
      w_err = U * w.abs
      w_up = w + w_err

      lo_f = lo.to_f
      hi_f = hi.to_f
      lo_conv = lo.is_a?(Float) ? 0.0 : U * lo_f.abs
      hi_conv = hi.is_a?(Float) ? 0.0 : U * hi_f.abs

      lo_term = lo_f * w
      hi_term = hi_f * w
      lo_err = (lo_conv * w_up) + (lo_f.abs * w_err) + (U * lo_term.abs)
      hi_err = (hi_conv * w_up) + (hi_f.abs * w_err) + (U * hi_term.abs)

      [lo_term, lo_err, hi_term, hi_err]
    end

    def add(lo_term, lo_err, hi_term, hi_err)
      accumulate_lo(lo_term, lo_err)
      accumulate_hi(hi_term, hi_err)
    end

    def remove(lo_term, lo_err, hi_term, hi_err)
      accumulate_lo(-lo_term, lo_err)
      accumulate_hi(-hi_term, hi_err)
    end

    def lower
      @sum_lo - @err_lo
    end

    def upper
      @sum_hi + @err_hi
    end

    def width
      upper - lower
    end

    private

    def accumulate_lo(term, err)
      @sum_lo += term
      @err_lo += err + (U * @sum_lo.abs)
    end

    def accumulate_hi(term, err)
      @sum_hi += term
      @err_hi += err + (U * @sum_hi.abs)
    end
  end
end
