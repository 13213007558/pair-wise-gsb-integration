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
    it 'should match per-component scalar integration with fewer callback evaluations' do
      components = (1..12).map { |k| lambda { |x| Math.exp(-k * x / 5.0) } }
      vector_evals = 0
      result = Integration.integrate_vector(0, Math::PI, tolerance: 1e-9) do |x|
        vector_evals += 1
        components.map { |c| c[x] }
      end

      scalar_evals = 0
      expected = components.map do |c|
        Integration.integrate(0, Math::PI, tolerance: 1e-9, method: :adaptive_quadrature) do |x|
          scalar_evals += 1
          c[x]
        end
      end

      result.values.each_with_index do |value, i|
        value.should be_within(1e-7).of(expected[i])
        value.should be_within(1e-7).of(5.0 * (1 - Math.exp(-(i + 1) * Math::PI / 5.0)) / (i + 1))
      end
      result.evaluations.should == vector_evals
      result.evaluations.should be < scalar_evals
    end

    it 'should only permute the results when components are permuted' do
      f = lambda { |x| [x**2, Math.sin(3 * x), Math.exp(-x), 1.0 / (1 + x)] }
      perm = [2, 0, 3, 1]
      options = { abs_tolerance: [1e-9, 1e-8, 1e-10, 1e-7], rel_tolerance: 1e-9 }

      original = Integration.integrate_vector(0, 2, options, &f)
      permuted = Integration.integrate_vector(0, 2, options.merge(
        abs_tolerance: perm.map { |i| options[:abs_tolerance][i] }
      )) { |x| perm.map { |i| f[x][i] } }

      perm.each_with_index do |i, j|
        permuted.values[j].should == original.values[i]
        permuted.diagnostics[j][:error].should == original.diagnostics[i][:error]
      end
      permuted.evaluations.should == original.evaluations
    end

    it 'should not let a large converged component mask a small spiky one' do
      spike = lambda do |x|
        [1.0, 1e-6 * Math.exp(-((x - 0.5) / 1e-3)**2)]
      end
      result = Integration.integrate_vector(0, 1,
        abs_tolerance: [1e-8, 1e-14], rel_tolerance: [1e-10, 1e-8], &spike)

      result.values[0].should be_within(1e-8).of(1.0)
      result.values[1].should be_within(1e-12).of(1e-6 * 1e-3 * Math.sqrt(Math::PI))
      result.diagnostics[1][:converged].should == true
      result.diagnostics[1][:error].should be <= 1e-12
    end

    it 'should return an empty result deterministically for an empty vector' do
      calls = 0
      f = lambda { |_x| calls += 1; [] }
      first = Integration.integrate_vector(0, 1, &f)
      second = Integration.integrate_vector(0, 1, &f)

      first.values.should == []
      first.diagnostics.should == []
      first.evaluations.should == 1
      first.converged?.should == true
      second.values.should == first.values
      second.evaluations.should == first.evaluations
    end

    it 'should raise when the callback changes length during integration' do
      calls = 0
      f = lambda do |x|
        calls += 1
        calls == 1 ? [x] : [x, x]
      end
      lambda { Integration.integrate_vector(0, 1, &f) }.should raise_exception(/changed length/)
    end

    it 'should raise when the callback does not return an Array' do
      lambda { Integration.integrate_vector(0, 1) { |x| x } }.should raise_exception(/Array/)
    end

    it 'should raise when a tolerance array does not match the vector length' do
      lambda do
        Integration.integrate_vector(0, 1, abs_tolerance: [1e-8, 1e-8]) { |x| [x] }
      end.should raise_exception(/abs_tolerance/)
    end
  end
end
