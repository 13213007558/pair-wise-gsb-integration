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
    unless Integration.has_gsl?
      skip 'rb-gsl is not installed'
    end
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qng }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
    Integration.integrate(0, 1, { tolerance: 1e-12, method: :qag }, &normal_pdf).should be_within(1e-11).of(0.341344746068)
  end

  it 'keeps the finite integration interface numeric' do
    result = Integration.integrate(1, 2, method: :simpson, tolerance: 1e-8) { |x| x**2 }
    result.should be_an(Float)
    result.should be_within(1e-6).of(7.0 / 3)
  end

  it 'proves an exponentially decaying tail with a tail bound' do
    f = lambda { |x| Math.exp(-x) }
    result = Integration.integrate_infinite(
      0,
      Integration::Infinity,
      tolerance: 1e-10,
      tail_bound: ->(x) { x.zero? ? 1.0 : Math.exp(-x) },
      &f
    )

    result.should be_an(Integration::InfiniteIntegralResult)
    result.should be_within(1e-9).of(1)
    result.interval_error.should be_finite
    result.tail_bound.should be <= 1e-10
    result.unproven.should == 0
    result.should be_proven
  end

  it 'proves a slowly decaying integrable power tail' do
    f = lambda { |x| 1.0 / (x * Math.sqrt(x)) }
    result = Integration.integrate(
      1,
      Integration::Infinity,
      tolerance: 1e-8,
      tail_bound: ->(x) { 2.0 / Math.sqrt(x) },
      &f
    )

    result.should be_within(result.tail_bound + 1e-7).of(2)
    result.interval_error.should be_finite
    result.tail_bound.should be_finite
    result.tail_bound.should be <= 1e-8
    result.unproven.should == 0
    result.should be_proven
    result.iterations.should > 5
    result.upper_cutoff.should > 1e10
  end

  it 'reports an unproven divergent tail instead of choosing a large endpoint' do
    result = Integration.integrate(
      1,
      Integration::Infinity,
      cutoff_factor: 2,
      max_iterations: 8
    ) { |x| 1.0 / x }

    result.should be_within(1e-12).of(Math.log(2))
    result.interval_error.should be_finite
    result.tail_bound.should == Float::INFINITY
    result.unproven.should == Float::INFINITY
    result.should_not be_proven
    result.upper_cutoff.should == 2
  end

  it 'supports an infinite lower bound with a separate tail bound' do
    f = lambda { |x| -1.0 / (x.abs * Math.sqrt(x.abs)) }
    result = Integration.integrate(
      Integration::MInfinity,
      -1,
      tolerance: 1e-8,
      lower_tail_bound: ->(x) { 2.0 / Math.sqrt(x.abs) },
      &f
    )

    result.should be_within(1e-7).of(-2)
    result.unproven.should == 0
  end

  it 'supports doubly infinite bounds with a direction-aware tail bound' do
    normal_pdf = lambda { |x| (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-(x**2 / 2)) }
    tail_bound = lambda do |x, direction|
      expected_sign = direction == :lower ? -1 : 1
      x.should satisfy { |value| (value <=> 0) == expected_sign }
      normal_pdf.call(x) / x.abs
    end
    result = Integration.integrate(
      Integration::MInfinity,
      Integration::Infinity,
      tolerance: 1e-10,
      tail_bound: tail_bound,
      &normal_pdf
    )

    result.should be_within(1e-9).of(1)
    result.tail_bound.should be_finite
    result.unproven.should == 0
  end

  it 'reports a finite computed part when no tail bound is supplied' do
    result = Integration.integrate_infinite(
      0,
      Integration::Infinity,
      cutoff: 2
    ) { |x| Math.exp(-x) }

    result.should be_within(1e-10).of(1 - Math.exp(-2))
    result.interval_error.should be_finite
    result.tail_bound.should == Float::INFINITY
    result.unproven.should == Float::INFINITY
    result.should_not be_proven
  end

  it 'supports a custom initial cutoff with a tail bound' do
    result = Integration.integrate_infinite(
      1,
      Integration::Infinity,
      cutoff: 3,
      tolerance: 1e-8,
      tail_bound: ->(x) { 2.0 / Math.sqrt(x) }
    ) { |x| 1.0 / (x * Math.sqrt(x)) }

    result.should be_within(result.tail_bound + 1e-7).of(2)
    result.upper_cutoff.should >= 3
    result.unproven.should == 0
  end

  it 'rejects a negative tail bound' do
    expect do
      Integration.integrate_infinite(
        0,
        Integration::Infinity,
        tail_bound: ->(_x) { -1 }
      ) { |x| Math.exp(-x) }
    end.to raise_error(RuntimeError, /nonnegative number/)
  end
end
