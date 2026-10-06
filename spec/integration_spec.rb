$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe Integration do
  a = lambda { |x| x**2 }
  b = lambda { |x| Math.log(x) / x**2 }
  b2 = lambda { |x| -(Math.log(x) + 1) / x }

  # Integration over [1,2]=x^3/3=7/3
  methods = [:rectangle, :trapezoid, :simpson, :adaptive_quadrature,  :romberg,
             :gauss, :gauss_kronrod, :simpson3by8, :boole, :open_trapezoid,
             :milne]

  methods.each do |m|
    it "should integrate int_{1}^2{2} x^2 correctly with ruby method #{m}" do
      Integration.integrate(1, 2, { method: m, tolerance: 1e-8 }, &a).should be_within(1e-6).of(7.0 / 3)
    end

    it "should integrate int_{1}^2{2} log(x)/x^2 correctly with ruby method #{m}" do
      Integration.integrate(1, 2, { method: m, tolerance: 1e-8 }, &b).should be_within(1e-6).of(
        b2[2] - b2[1]
      )
    end
  end

  it 'should return correct for trapezoid' do
    a = rand
    b = rand * 10
    f = lambda { |x| x * a }
    Integration.trapezoid(0, b, 2, &f).should be_within(1e-14).of((a * b**2) / 2.0)
  end

  it 'should return a correct value for a complex integration with ruby methods' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :simpson }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :adaptive_quadrature }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
  end

  it 'should return a correct value for a complex integration with gsl methods' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qng }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qag }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
  end

  it 'should return correct integration for infinity bounds' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }

    Integration.integrate(Integration::MInfinity, Integration::Infinity, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(1)
  end

  it 'should return correct integration for infinity lower bound' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(Integration::MInfinity, 0, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(0.5)
  end

  it 'should return correct integration for infinity upper bound' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, Integration::Infinity, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(0.5)
  end

  it 'should raise an error if a ruby methods is called with infinite bounds' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    lambda { Integration.integrate(0, Integration::Infinity, { method: :simpson }, &normal_pdf).should be_within(1e-09).of(0.5) }.should raise_exception
  end

  describe '.validated_integrate' do
    let(:squared_enclosure) do
      lambda do |lower, upper|
        lower = Rational(lower)
        upper = Rational(upper)
        [lower**2, upper**2]
      end
    end

    let(:analytic_integral) do
      lambda { |x| Rational(x)**3 / 3 }
    end

    let(:centered_squared_enclosure) do
      lambda do |lower, upper|
        lower = Rational(lower)
        upper = Rational(upper)
        midpoint = (lower + upper) / 2
        half_width = (upper - lower) / 2
        slack = midpoint.abs * half_width + half_width**2 / 4
        [midpoint**2 - slack, midpoint**2 + slack]
      end
    end

    it 'keeps analytical truth inside increasingly tight Rational enclosures' do
      truth = analytic_integral.call(2) - analytic_integral.call(1)
      results = [1, 3, 7, 15, 31].map do |budget|
        Integration.validated_integrate(
          1r, 2r, tolerance: 0, max_evaluations: budget, &squared_enclosure
        )
      end

      results.each_cons(2) do |previous, current|
        expect(current.lower_bound).to be >= previous.lower_bound
        expect(current.upper_bound).to be <= previous.upper_bound
      end

      results.each do |result|
        expect(truth).to be_between(result.lower_bound, result.upper_bound)
        expect(result.state).to eq(:budget_exhausted)
      end
      expect(results.last.evaluations).to eq(31)
    end

    it 'returns the current enclosure when the evaluation budget is exhausted' do
      result = Integration.validated_integrate(
        1r, 2r, tolerance: 0, max_evaluations: 3, &squared_enclosure
      )

      expect(result.state).to eq(:budget_exhausted)
      expect(result.converged?).to be(false)
      expect(result.evaluations).to eq(3)
      expect(result.interval_count).to eq(2)
      expect(result.lower_bound).to eq(13r / 8)
      expect(result.upper_bound).to eq(25r / 8)
      expect(result.error).to eq(3r / 2)
    end

    it 'converges with Rational bounds without sampling a scalar callback' do
      result = Integration.validated_integrate(
        1r, 2r, tolerance: 1r / 10, max_evaluations: 127,
        &centered_squared_enclosure
      )
      truth = analytic_integral.call(2) - analytic_integral.call(1)

      expect(result.state).to eq(:converged)
      expect(result.converged?).to be(true)
      expect(result.evaluations).to be <= 127
      expect(result.error).to be <= 1r / 10
      expect(truth).to be_between(result.lower_bound, result.upper_bound)
    end

    it 'requires explicit outward rounding for Float bounds' do
      expect do
        Integration.validated_integrate(1.0, 2.0, &squared_enclosure)
      end.to raise_error(ArgumentError, /explicit outward/)
    end

    it 'validates Float-bounded integration through explicit outward rounding' do
      result = Integration.validated_integrate(
        1.0, 2.0, tolerance: 0.1, max_evaluations: 127,
        outward_round: :float, &centered_squared_enclosure
      )
      truth = analytic_integral.call(2) - analytic_integral.call(1)

      expect(result.state).to eq(:converged)
      expect(result.lower_bound).to be_a(Float)
      expect(result.upper_bound).to be_a(Float)
      expect(result.lower_bound).to be <= truth
      expect(result.upper_bound).to be >= truth
      expect(result.error).to be <= 0.1
    end

    it 'accepts a callable outward rounding object' do
      rounder = lambda do |value, direction|
        nearest = value.to_f
        direction == :down ? nearest.prev_float : nearest.next_float
      end
      result = Integration.validated_integrate(
        1.0, 2.0, tolerance: 0.1, max_evaluations: 127,
        outward_round: rounder, &centered_squared_enclosure
      )

      expect(result.state).to eq(:converged)
      expect(7r / 3).to be_between(result.lower_bound, result.upper_bound)
    end

    it 'accepts a Range enclosure and preserves reversed-bound orientation' do
      oracle = lambda do |lower, upper|
        (Rational(lower)**2..Rational(upper)**2)
      end

      result = Integration.validated_integrate(
        2r, 1r, tolerance: 0, max_evaluations: 3, &oracle
      )

      expect(result.lower_bound).to eq(-25r / 8)
      expect(result.upper_bound).to eq(-13r / 8)
      expect(result.bounds).to eq([-25r / 8, -13r / 8])
      expect(result[:state]).to eq(:budget_exhausted)
    end

    it 'rejects an ordinary scalar Float sampling callback' do
      expect do
        Integration.validated_integrate(1r, 2r) { |point| point.to_f }
      end.to raise_error(ArgumentError, /two-sided numeric enclosure/)
    end

    it 'rejects non-finite oracle bounds' do
      expect do
        Integration.validated_integrate(1r, 2r) do
          [0.0, Float::NAN]
        end
      end.to raise_error(ArgumentError, /two-sided numeric enclosure/)
    end
  end
end
