$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe Integration do
  describe '.tanh_sinh' do
    it 'integrates log(x) over [0,1] with an integrable endpoint singularity' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| Math.log(x) })
      result.value.should be_within(1e-9).of(-1.0)
      result.should be_converged
      result.stop_reason.should eq :converged
    end

    it 'integrates 1/sqrt(x) over [0,1] with an integrable endpoint singularity' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| 1.0 / Math.sqrt(x) })
      result.value.should be_within(1e-9).of(2.0)
      result.should be_converged
    end

    it 'integrates smooth functions over finite intervals' do
      result = Integration.tanh_sinh(1, 2, &lambda { |x| x**2 })
      result.value.should be_within(1e-9).of(7.0 / 3)
      result.should be_converged
    end

    it 'handles reversed bounds' do
      result = Integration.tanh_sinh(1, 0, &lambda { |x| Math.log(x) })
      result.value.should be_within(1e-9).of(1.0)
    end

    it 'never evaluates the integrand at the endpoints' do
      endpoints = [0.0, 1.0]
      f = lambda do |x|
        fail "evaluated at endpoint #{x}" if endpoints.include? x
        Math.log(x) / Math.sqrt(x)
      end
      result = Integration.tanh_sinh(0, 1, &f)
      result.value.should be_within(1e-8).of(-4.0)
      result.should be_converged
    end

    it 'reports the number of evaluations and the stop reason' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| Math.log(x) })
      result.evaluations.should be > 0
      result.levels.should be >= 1
      [:converged, :max_levels, :max_evaluations, :non_finite].should include result.stop_reason
    end

    it 'checks discretization error and tail truncation separately' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| Math.log(x) })
      result.discretization_error.should be <= 1e-10
      result.tail_error.should be <= 1e-10
      result.error.should eq [result.discretization_error, result.tail_error].max
    end

    it 'does not claim convergence for the divergent integral of 1/x' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| 1.0 / x })
      result.should_not be_converged
      result.stop_reason.should_not eq :converged
    end

    it 'does not claim convergence for the divergent integral of 1/x^2' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| 1.0 / x**2 })
      result.should_not be_converged
      result.stop_reason.should_not eq :converged
    end

    it 'stops with :max_evaluations when the evaluation budget is exhausted' do
      result = Integration.tanh_sinh(0, 1, max_evaluations: 10, &lambda { |x| Math.log(x) })
      result.stop_reason.should eq :max_evaluations
      result.should_not be_converged
      result.evaluations.should be <= 10
    end

    it 'stops with :max_levels when the level limit is reached first' do
      result = Integration.tanh_sinh(0, 1, max_levels: 1, tolerance: 1e-30,
                                     &lambda { |x| Math.log(x) })
      result.stop_reason.should eq :max_levels
      result.should_not be_converged
    end

    it 'stops with :non_finite when the integrand yields NaN' do
      result = Integration.tanh_sinh(0, 1, &lambda { |x| x > 0.5 ? Float::NAN : 1.0 })
      result.stop_reason.should eq :non_finite
      result.should_not be_converged
    end

    it 'returns a zero result for an empty interval' do
      result = Integration.tanh_sinh(1, 1, &lambda { |x| x })
      result.value.should eq 0.0
      result.should be_converged
      result.evaluations.should eq 0
    end

    it 'rejects infinite bounds' do
      lambda do
        Integration.tanh_sinh(0, Integration::Infinity, &lambda { |x| x })
      end.should raise_exception(/infinite|Non-numeric/)
      lambda do
        Integration.tanh_sinh(0, Float::INFINITY, &lambda { |x| x })
      end.should raise_exception(/infinite/)
    end
  end

  describe '.integrate with method: :tanh_sinh' do
    it 'integrates log(x) over [0,1] when explicitly selected' do
      value = Integration.integrate(0, 1, method: :tanh_sinh, &lambda { |x| Math.log(x) })
      value.should be_within(1e-9).of(-1.0)
    end

    it 'integrates 1/sqrt(x) over [0,1] when explicitly selected' do
      value = Integration.integrate(0, 1, method: :tanh_sinh, &lambda { |x| 1.0 / Math.sqrt(x) })
      value.should be_within(1e-9).of(2.0)
    end

    it 'raises instead of returning a value for a divergent integral' do
      lambda do
        Integration.integrate(0, 1, method: :tanh_sinh, &lambda { |x| 1.0 / x })
      end.should raise_exception(/did not converge/)
    end

    it 'keeps the default method unchanged when no method is given' do
      if Integration.has_gsl?
        expect(Integration).to receive(:integrate_gsl).with(1.0, 2.0, hash_including(method: :qag)).and_call_original
      else
        expect(Integration).to receive(:integrate_ruby).with(1.0, 2.0, hash_including(method: :simpson)).and_call_original
      end
      Integration.integrate(1, 2, &lambda { |x| x**2 }).should be_within(1e-6).of(7.0 / 3)
    end
  end
end
