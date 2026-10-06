$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe Integration do
  describe '.validated_integrate' do
    # Natural, exactly computable enclosure of x^2 on [0, 1]:
    # x^2 is increasing there, so [a^2, b^2] bounds it on [a, b].
    let(:square_enclosure) { lambda { |a, b| [a * a, b * b] } }

    # Wide but valid enclosure of x^2 on [a, b] subset of [0, 1]:
    # centered at the midpoint value with half-width 2 * (b - a).
    # It tightens monotonically as subintervals shrink.
    let(:wide_square_oracle) do
      lambda do |a, b|
        center = ((a + b) / 2)**2
        half_width = 2 * (b - a)
        [center - half_width, center + half_width]
      end
    end

    it 'encloses the true value of int x^2 with exact Rational bounds' do
      result = Integration.validated_integrate(Rational(0), Rational(1),
                                               tolerance: Rational(1, 1000),
                                               &square_enclosure)

      expect(result.status).to eq(:converged)
      expect(result).to be_converged
      expect(result.lower).to be_a(Rational)
      expect(result.upper).to be_a(Rational)
      expect(result.lower).to be <= Rational(1, 3)
      expect(result.upper).to be >= Rational(1, 3)
      expect(result.width).to be <= Rational(1, 1000)
      expect(result).to include(Rational(1, 3))
    end

    it 'is exact for a constant function with Rational enclosure' do
      result = Integration.validated_integrate(Rational(0), Rational(2),
                                               tolerance: Rational(1, 10**6)) { |_a, _b| [1, 1] }

      expect(result.status).to eq(:converged)
      expect(result.lower).to eq(Rational(2))
      expect(result.upper).to eq(Rational(2))
      expect(result.width).to eq(0)
    end

    it 'tightens monotonically as the evaluation budget grows' do
      budgets = [5, 21, 101, 501]
      results = budgets.map do |budget|
        Integration.validated_integrate(Rational(0), Rational(1),
                                        tolerance: 0,
                                        max_evaluations: budget,
                                        &wide_square_oracle)
      end

      results.each do |result|
        expect(result.status).to eq(:budget_exhausted)
        expect(result).to include(Rational(1, 3))
      end

      lowers = results.map(&:lower)
      uppers = results.map(&:upper)
      widths = results.map(&:width)

      expect(lowers.each_cons(2).all? { |a, b| b >= a }).to be(true)
      expect(uppers.each_cons(2).all? { |a, b| b <= a }).to be(true)
      expect(widths.each_cons(2).all? { |a, b| b <= a }).to be(true)
    end

    it 'keeps the true value enclosed at every refinement step' do
      (1..8).each do |splits|
        budget = (2 * splits) + 1
        result = Integration.validated_integrate(Rational(0), Rational(1),
                                                 tolerance: 0,
                                                 max_evaluations: budget,
                                                 &wide_square_oracle)
        expect(result).to include(Rational(1, 3))
      end
    end

    it 'encloses e - 1 for exp with an outward-rounded Float enclosure' do
      # Explicit outward rounding on the caller side: exp is increasing,
      # and the relative slack covers the rounding of Math.exp.
      enclosure = lambda do |a, b|
        [Math.exp(a) * (1 - 1e-13), Math.exp(b) * (1 + 1e-13)]
      end

      result = Integration.validated_integrate(0.0, 1.0,
                                               tolerance: 1e-3,
                                               &enclosure)

      expect(result.status).to eq(:converged)
      expect(result.lower).to be <= Math::E - 1
      expect(result.upper).to be >= Math::E - 1
      expect(result.width).to be <= 1e-3
    end

    it 'encloses pi/4 for 1/(1+x^2) with exact Rational enclosure' do
      # 1/(1+x^2) is decreasing on [0, 1].
      enclosure = lambda do |a, b|
        [Rational(1, 1 + (b * b)), Rational(1, 1 + (a * a))]
      end

      result = Integration.validated_integrate(Rational(0), Rational(1),
                                               tolerance: Rational(1, 500),
                                               &enclosure)

      expect(result.status).to eq(:converged)
      expect(result.lower).to be_a(Rational)
      expect(result.lower.to_f).to be <= Math::PI / 4
      expect(result.upper.to_f).to be >= Math::PI / 4
    end

    it 'returns the current enclosure and :budget_exhausted when the budget runs out' do
      result = Integration.validated_integrate(Rational(0), Rational(1),
                                               tolerance: Rational(1, 10**9),
                                               max_evaluations: 11,
                                               &square_enclosure)

      expect(result.status).to eq(:budget_exhausted)
      expect(result).not_to be_converged
      expect(result.evaluations).to be <= 11
      expect(result.width).to be > Rational(1, 10**9)
      expect(result).to include(Rational(1, 3))
    end

    it 'never exceeds the evaluation budget' do
      result = Integration.validated_integrate(0.0, 1.0,
                                               tolerance: 1e-12,
                                               max_evaluations: 101,
                                               &wide_square_oracle)

      expect(result.evaluations).to be <= 101
      expect(result.status).to eq(:budget_exhausted)
      expect(result).to include(Rational(1, 3))
    end

    it 'rounds Float results outward so the true value stays inside' do
      # Constant function 1 on [0, 1]: exact integral is 1. The Float
      # accumulation must still bracket it despite rounding.
      result = Integration.validated_integrate(0.0, 1.0,
                                               tolerance: 1e-12,
                                               max_evaluations: 2001) { |_a, _b| [1.0, 1.0] }

      expect(result.lower).to be <= 1.0
      expect(result.upper).to be >= 1.0
      expect(result.width).to be < 1e-9
    end

    it 'refuses Float enclosures in exact arithmetic mode' do
      floaty = lambda { |a, b| [Math.exp(a.to_f), Math.exp(b.to_f)] }

      expect do
        Integration.validated_integrate(Rational(0), Rational(1),
                                        arithmetic: :exact,
                                        &floaty)
      end.to raise_error(/exact/)
    end

    it 'falls back to outward-rounded Float arithmetic for Float enclosures' do
      floaty = lambda { |a, b| [Math.exp(a) * (1 - 1e-13), Math.exp(b) * (1 + 1e-13)] }

      result = Integration.validated_integrate(0, 1, tolerance: 1e-3, &floaty)

      expect(result.status).to eq(:converged)
      expect(result).to include(Math::E - 1)
    end

    it 'rejects enclosures whose lower bound exceeds the upper bound' do
      bad = lambda { |_a, _b| [1, 0] }

      expect do
        Integration.validated_integrate(0, 1, &bad)
      end.to raise_error(/Invalid enclosure/)
    end

    it 'rejects enclosures that do not return a pair' do
      bad = lambda { |_a, _b| 42 }

      expect do
        Integration.validated_integrate(0, 1, &bad)
      end.to raise_error(/pair/)
    end

    it 'requires an enclosure block' do
      expect do
        Integration.validated_integrate(0, 1)
      end.to raise_error(/No enclosure/)
    end

    it 'rejects infinite bounds' do
      expect do
        Integration.validated_integrate(0, Integration::Infinity) { |a, b| [a, b] }
      end.to raise_error(/infinite/)
    end

    it 'rejects non-numeric bounds' do
      expect do
        Integration.validated_integrate('a', 1) { |a, b| [a, b] }
      end.to raise_error(/Non-numeric/)
    end

    it 'exposes midpoint and to_a helpers' do
      result = Integration.validated_integrate(Rational(0), Rational(1),
                                               tolerance: Rational(1, 100),
                                               &square_enclosure)

      expect(result.to_a).to eq([result.lower, result.upper])
      expect(result.midpoint).to eq((result.lower + result.upper) / 2)
      expect(result.evaluations).to be_an(Integer)
    end
  end
end
