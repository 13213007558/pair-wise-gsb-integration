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

  describe '.integrate_vector' do
    let(:amplitude) { 1e-3 }
    let(:width) { 0.002 }
    let(:center) { 0.75 }
    let(:component_count) { 12 }

    let(:spike) do
      lambda do |x|
        scaled = (x - center) / width
        amplitude * Math.exp(-(scaled * scaled))
      end
    end

    let(:large_component) { lambda { |_x| 1000.0 } }
    let(:callback) do
      lambda do |x|
        Array.new(component_count - 1) { large_component.call(x) } + [spike.call(x)]
      end
    end
    let(:abs_tolerances) { Array.new(component_count - 1, 1e-6) + [1e-10] }
    let(:spike_integral) { amplitude * width * Math.sqrt(Math::PI) }

    it 'shares samples and keeps per-component diagnostics independent' do
      callback_count = 0
      counted_callback = lambda do |x|
        callback_count += 1
        callback.call(x)
      end

      result = Integration.integrate_vector(
        0, 1,
        abs_tolerance: abs_tolerances,
        relative_tolerance: 0,
        &counted_callback
      )

      expect(result).to be_converged
      expect(result.values.first).to be_within(1e-12).of(1000.0)
      expect(result.values.last).to be_within(1e-10).of(spike_integral)
      expect(result.errors.last).to be <= 1e-10
      expect(result.callback_count).to eq(callback_count)
      expect(result.components.length).to eq(component_count)
      expect(result.components.map(&:converged)).to all(eq(true))
      expect(result.components.last.abs_tolerance).to eq(1e-10)
    end

    it 'uses fewer bundled callbacks than integrating every component separately' do
      bundled_count = 0
      bundled_result = Integration.integrate_vector(
        0, 1,
        abs_tolerance: abs_tolerances,
        relative_tolerance: 0
      ) do |x|
        bundled_count += 1
        callback.call(x)
      end

      separate_counts = Array.new(component_count, 0)
      separate_results = Array.new(component_count)
      (0...(component_count - 1)).each do |index|
        separate_results[index] = Integration.adaptive_quadrature(0, 1, 1e-6) do |x|
          separate_counts[index] += 1
          large_component.call(x)
        end
      end
      separate_results[-1] = Integration.adaptive_quadrature(0, 1, 1e-10) do |x|
        separate_counts[-1] += 1
        spike.call(x)
      end

      expect(bundled_count).to be < separate_counts.reduce(:+)
      bundled_result.values.each_with_index do |value, index|
        tolerance = abs_tolerances[index] * 10
        expect(value).to be_within(tolerance).of(separate_results[index])
      end
    end

    it 'reorders values, tolerances, and diagnostics with component order' do
      order = (0...component_count).to_a.reverse
      permuted_callback = lambda do |x|
        values = callback.call(x)
        order.map { |index| values[index] }
      end
      permuted_tolerances = order.map { |index| abs_tolerances[index] }

      original = Integration.integrate_vector(0, 1, abs_tolerance: abs_tolerances,
                                                  relative_tolerance: 0, &callback)
      permuted = Integration.integrate_vector(0, 1, abs_tolerance: permuted_tolerances,
                                                    relative_tolerance: 0,
                                                    &permuted_callback)

      expect(permuted.values).to eq(order.map { |index| original.values[index] })
      expect(permuted.errors).to eq(order.map { |index| original.errors[index] })
      expect(permuted.components.map(&:index)).to eq((0...component_count).to_a)
      expect(permuted.callback_count).to eq(original.callback_count)
      permuted.components.each_with_index do |component, index|
        original_component = original.components[order[index]]
        expect(component.value).to eq(original_component.value)
        expect(component.error).to eq(original_component.error)
        expect(component.abs_tolerance).to eq(original_component.abs_tolerance)
        expect(component.relative_tolerance).to eq(original_component.relative_tolerance)
        expect(component.target_error).to eq(original_component.target_error)
        expect(component.converged?).to eq(original_component.converged?)
      end
    end

    it 'returns an empty result for an empty vector' do
      calls = 0
      result = Integration.integrate_vector(0, 1, abs_tolerance: []) do |_x|
        calls += 1
        []
      end

      expect(result).to be_converged
      expect(result.values).to eq([])
      expect(result.errors).to eq([])
      expect(result.components).to eq([])
      expect(result.callback_count).to eq(calls)
      expect(result.callback_count).to eq(1)
    end

    it 'raises before returning when callback length changes' do
      calls = 0
      expect do
        Integration.integrate_vector(0, 1) do |x|
          calls += 1
          x == 0 ? [1.0] : [1.0, 2.0]
        end
      end.to raise_error(ArgumentError, /changed length from 1 to 2/)
      expect(calls).to be > 1
    end

    it 'preserves the sign of reversed bounds' do
      forward = Integration.integrate_vector(0, 1, abs_tolerance: 1e-10) { |x| [x, x**2] }
      backward = Integration.integrate_vector(1, 0, abs_tolerance: 1e-10) { |x| [x, x**2] }

      expect(backward.values).to eq(forward.values.map { |value| -value })
      expect(backward.errors).to eq(forward.errors)
    end

    it 'uses separate absolute and relative targets for components' do
      result = Integration.integrate_vector(
        0, 1,
        abs_tolerance: [1e-8, 1e-12],
        relative_tolerance: [1e-10, 1e-9]
      ) { |x| [x, x**2] }

      expect(result.components[0].target_error).to be_within(1e-16).of(1e-8 + 1e-10 * 0.5)
      expect(result.components[1].target_error).to be_within(1e-16).of(1e-12 + 1e-9 / 3.0)
      expect(result).to be_converged
    end
  end
end
