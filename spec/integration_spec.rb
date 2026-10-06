$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe Integration do
  def complex_monomial_moment(degree, omega)
    endpoint = Complex(Math.cos(omega), Math.sin(omega))
    return (endpoint - 1.0) / Complex(0.0, omega) if degree.zero?

    endpoint / Complex(0.0, omega) - degree.to_f / Complex(0.0, omega) *
      complex_monomial_moment(degree - 1, omega)
  end

  def polynomial_oscillatory_integral(omega, oscillation)
    moment = complex_monomial_moment(0, omega) +
             2.0 * complex_monomial_moment(1, omega) +
             3.0 * complex_monomial_moment(2, omega)
    oscillation == :sin ? moment.imag : moment.real
  end

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

  it 'should integrate a polynomial amplitude with cosine using analytically exact moments' do
    omega = 1000.0
    evaluations = 0
    result = Integration.oscillatory(0, 1, omega, oscillation: :cos) do |x|
      evaluations += 1
      1 + (2 * x) + (3 * x * x)
    end
    expected = polynomial_oscillatory_integral(omega, :cos)

    expect(result).to be_a(Integration::OscillatoryResult)
    expect(result.value).to be_within(1e-12).of(expected)
    expect(result.error).to be <= 1e-8
    expect(result.amplitude_evaluations).to eq(evaluations)
    expect(result.converged).to be(true)
  end

  it 'should integrate a polynomial amplitude with sine using analytically exact moments' do
    omega = 1000.0
    result = Integration.oscillatory(0, 1, omega, oscillation: :sin) do |x|
      1 + (2 * x) + (3 * x * x)
    end
    expected = polynomial_oscillatory_integral(omega, :sin)

    expect(result.value).to be_within(1e-12).of(expected)
    expect(result.error).to be <= 1e-8
  end

  it 'should not let amplitude evaluations grow linearly with frequency' do
    frequencies = [10.0, 100.0, 1000.0, 10000.0, 100000.0]
    evaluation_counts = frequencies.map do |omega|
      count = 0
      result = Integration.oscillatory(0, 1, omega, oscillation: :cos) do |x|
        count += 1
        1 + (2 * x) + (3 * x * x)
      end
      expected = polynomial_oscillatory_integral(omega, :cos)

      expect(result.value).to be_within(1e-10).of(expected)
      count
    end

    expect(evaluation_counts.max).to eq(evaluation_counts.min)
    expect(evaluation_counts.last).to be < 32
  end

  it 'should fall back to ordinary quadrature when the frequency is near zero' do
    evaluations = 0
    result = Integration.oscillatory(0, 1, 0.0, tolerance: 1e-10) do |x|
      evaluations += 1
      1 + (2 * x) + (3 * x * x)
    end

    expect(result.value).to be_within(1e-10).of(3)
    expect(result.error).to be <= 1e-10
    expect(result.amplitude_evaluations).to eq(evaluations)
  end

  it 'should handle a small sine frequency without cancellation' do
    result = Integration.oscillatory(0, 1, 1e-9, oscillation: :sin) do |x|
      1 + (2 * x) + (3 * x * x)
    end

    expect(result.value).to be_within(1e-12).of(23.0 / 12 * 1e-9)
  end

  it 'should return zero for a zero-width oscillatory interval without evaluating the amplitude' do
    evaluations = 0
    result = Integration.oscillatory(2, 2, 1000.0) do
      evaluations += 1
      1
    end

    expect(result.value).to eq(0)
    expect(result.error).to eq(0)
    expect(result.amplitude_evaluations).to eq(0)
    expect(evaluations).to eq(0)
  end

  it 'should reject non-finite floating point frequency' do
    expect do
      Integration.oscillatory(0, 1, Float::INFINITY) { |x| x }
    end.to raise_error(RuntimeError, 'Finite bounds and frequency are required')
  end

  it 'should preserve reversed finite bounds' do
    result = Integration.oscillatory(1, 0, 1000.0, oscillation: :sin) do |x|
      1 + (2 * x) + (3 * x * x)
    end
    expected = -polynomial_oscillatory_integral(1000.0, :sin)

    expect(result.value).to be_within(1e-12).of(expected)
  end

  it 'should reject infinite bounds for oscillatory integration' do
    expect do
      Integration.oscillatory(0, Integration::Infinity, 1.0) { |x| x }
    end.to raise_error(RuntimeError, /Finite bounds are required|Non-numeric bounds/)
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
end
