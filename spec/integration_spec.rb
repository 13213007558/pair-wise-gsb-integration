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
    next unless Integration.has_gsl?

    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qng }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qag }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
  end

  it 'should return correct integration for infinity bounds with GSL' do
    next unless Integration.has_gsl?

    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }

    Integration.integrate(Integration::MInfinity, Integration::Infinity, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(1)
  end

  it 'should return correct integration for infinity lower bound with GSL' do
    next unless Integration.has_gsl?

    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(Integration::MInfinity, 0, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(0.5)
  end

  it 'should return correct integration for infinity upper bound with GSL' do
    next unless Integration.has_gsl?

    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, Integration::Infinity, { tolerance: 1e-10 }, &normal_pdf).should be_within(1e-09).of(0.5)
  end

  it 'proves an exponential upper tail with a caller supplied tail bound' do
    f = lambda { |x| Math.exp(-x) }
    tail_bound = lambda { |cutoff| Math.exp(-cutoff) }
    result = Integration.integrate_infinite(
      0, Integration::Infinity,
      tolerance: 1e-8,
      tail_bound: tail_bound,
      &f
    )

    result.value.should be_within(1e-7).of(1)
    result.interval_error.should be <= 1e-8
    result.tail_bound.should be <= 1e-8
    result.error_bound.should be <= 1e-8
    result.unproven_part.should eq(0)
    result.should be_proven
  end

  it 'proves an integrable power tail despite slow decay' do
    f = lambda { |x| x**Rational(-3, 2) }
    tail_bound = lambda { |cutoff| 2 / Math.sqrt(cutoff) }
    result = Integration.integrate_infinite(
      1, Integration::Infinity,
      tolerance: 1e-5,
      tail_bound: tail_bound,
      max_extensions: 128,
      &f
    )

    result.value.should be_within(1e-4).of(2)
    result.interval_error.should be <= 1e-5
    result.tail_bound.should be <= 1e-5
    result.error_bound.should be <= 1e-5
    result.unproven_part.should eq(0)
    result.should be_proven
  end

  it 'reports an unproven infinite tail for a divergent comparison' do
    f = lambda { |x| 1 / x }
    finite_looking = Integration.adaptive_quadrature(1, 10_000, 1e-8) { |x| f.call(x) }
    result = Integration.integrate_infinite(
      1, Integration::Infinity,
      initial_cutoff: 10_000,
      tail_bound: ->(_cutoff) { nil },
      &f
    )

    finite_looking.should be_within(1e-6).of(Math.log(10_000))
    result.value.should be_within(1e-7).of(finite_looking)
    result.interval.should eq([1, 10_000])
    result.interval_error.should be <= 1e-8
    result.tail_bound.should eq(Float::INFINITY)
    result.unproven_part.should eq(Float::INFINITY)
    result.should_not be_proven
  end

  it 'does not mistake a stable finite sample for a convergent slow tail' do
    integrable_but_slow = lambda { |x| x**Rational(-3, 2) }
    finite_estimate = Integration.adaptive_quadrature(1, 10_000, 1e-8) { |x| integrable_but_slow.call(x) }
    result = Integration.integrate_infinite(
      1, Integration::Infinity,
      initial_cutoff: 10_000,
      &integrable_but_slow
    )

    finite_estimate.should be_within(1e-6).of(2 - 0.02)
    result.value.should be_within(1e-7).of(finite_estimate)
    result.tail_bound.should eq(Float::INFINITY)
    result.unproven_part.should eq(Float::INFINITY)
    result.should_not be_proven
  end

  it 'supports lower and upper infinite bounds with directional tail bounds' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2.0)) }
    bound = lambda do |cutoff, direction|
      x = cutoff.abs
      Math.exp(-(x**2 / 2.0)) / (x * Math.sqrt(2 * Math::PI))
    end
    options = {
      tolerance: 1e-7,
      split_at: 0,
      lower_tail_bound: bound,
      upper_tail_bound: bound
    }
    result = Integration.integrate_infinite(
      Integration::MInfinity, Integration::Infinity,
      options,
      &normal_pdf
    )

    result.value.should be_within(1e-6).of(1)
    result.interval.first.should be < 0
    result.interval.last.should be > 0
    result.error_bound.should be <= 1e-7
    result.unproven_part.should eq(0)
    result.should be_proven
  end

  it 'keeps the finite interval interface numeric' do
    Integration.integrate(1, 2, method: :simpson) { |x| x**2 }.should be_within(1e-8).of(7.0 / 3)
  end

  it 'requires a tail bound before returning a scalar infinite Ruby integral' do
    expect do
      Integration.integrate(0, Integration::Infinity, method: :simpson) { |x| Math.exp(-x) }
    end.to raise_error(RuntimeError, /tail bound/)
  end
end
