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

  def oscillatory_polynomial_integral(kind, frequency)
    w = frequency.to_f
    return 1.0 / 3.0 if w.zero? && [:cosine, :cos].include?(kind)
    return 0.0 if w.zero?

    if w.abs < 0.25
      if [:sine, :sin].include?(kind)
        (0..12).inject(0.0) do |sum, n|
          factorial = (1..(2 * n + 1)).inject(1, :*)
          sum + ((-1)**n) * w**(2 * n + 1) / (factorial * (2 * n + 4).to_f)
        end
      else
        (0..12).inject(0.0) do |sum, n|
          factorial = (1..(2 * n)).inject(1, :*)
          sum + ((-1)**n) * w**(2 * n) / (factorial * (2 * n + 3).to_f)
        end
      end
    elsif [:sine, :sin].include?(kind)
      -Math.cos(w) / w + 2 * Math.sin(w) / w**2 +
        2 * Math.cos(w) / w**3 - 2.0 / w**3
    else
      Math.sin(w) / w + 2 * Math.cos(w) / w**2 -
        2 * Math.sin(w) / w**3
    end
  end

  it 'provides an oscillatory quadrature result with diagnostics' do
    amplitude = lambda { |x| x**2 }
    result = Integration.oscillatory_integral(0, 1, 10, kind: :cosine, tolerance: 1e-9, &amplitude)

    expect(result.value).to be_within(1e-9).of(oscillatory_polynomial_integral(:cosine, 10))
    expect(result.error_estimate).to be >= 0
    expect(result.error_estimate).to be <= 1e-9
    expect(result.amplitude_evaluations).to be > 0
    expect(result.evaluations).to eq(result.amplitude_evaluations)
    expect(result.to_f).to eq(result.value)
  end

  it 'integrates polynomial amplitudes against both oscillatory factors' do
    amplitude = lambda { |x| x**2 }

    [:cosine, :sine].each do |kind|
      [2, 10, 100].each do |frequency|
        result = Integration.oscillatory_integral(
          0, 1, frequency, kind: kind, tolerance: 1e-10, &amplitude
        )

        expect(result.value).to be_within(1e-9).of(
          oscillatory_polynomial_integral(kind, frequency)
        )
      end
    end
  end

  it 'falls back stably when the frequency is near zero' do
    amplitude = lambda { |x| x**2 }

    cosine = Integration.oscillatory_integral(0, 1, 0, kind: :cosine, tolerance: 1e-10, &amplitude)
    sine = Integration.oscillatory_integral(0, 1, 0, kind: :sine, tolerance: 1e-10, &amplitude)
    tiny_sine = Integration.oscillatory_integral(
      0, 1, 1e-8, kind: :sine, tolerance: 1e-10, &amplitude
    )

    expect(cosine.value).to be_within(1e-12).of(1.0 / 3.0)
    expect(sine.value).to eq(0.0)
    expect(tiny_sine.value).to be_within(1e-12).of(2.5e-9)
    expect(cosine.error_estimate).to be >= 0
    expect(sine.error_estimate).to eq(0.0)
  end

  it 'does not let amplitude evaluations grow linearly with frequency' do
    amplitude = lambda { |x| x**2 }
    frequencies = [2, 1000, 100_000]

    counts = frequencies.map do |frequency|
      result = Integration.oscillatory_integral(
        0, 1, frequency, kind: :cosine, tolerance: 1e-9, &amplitude
      )

      expect(result.value).to be_within(1e-8).of(
        oscillatory_polynomial_integral(:cosine, frequency)
      )
      result.amplitude_evaluations
    end

    expect(counts.max).to be <= 2 * counts.min
    expect(counts.max).to be < 100
  end

  it 'reverses the integral value for reversed finite bounds' do
    amplitude = lambda { |x| x**2 }

    forward = Integration.oscillatory_integral(0, 1, 17, kind: :sine, &amplitude)
    backward = Integration.oscillatory_integral(1, 0, 17, kind: :sine, &amplitude)

    expect(backward.value).to be_within(1e-12).of(-forward.value)
  end
end
