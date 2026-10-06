$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe Integration do
  describe '.tanh_sinh' do
    it 'integrates log(x) over [0, 1] with an integrable endpoint singularity' do
      result = Integration.tanh_sinh(0, 1) { |x| Math.log(x) }
      expect(result).to be_converged
      expect(result.status).to eq(:ok)
      expect(result.value).to be_within(1e-9).of(-1.0)
    end

    it 'integrates 1/sqrt(x) over [0, 1] with an integrable endpoint singularity' do
      result = Integration.tanh_sinh(0, 1) { |x| 1.0 / Math.sqrt(x) }
      expect(result).to be_converged
      expect(result.value).to be_within(1e-9).of(2.0)
    end

    it 'never evaluates the integrand at the endpoints' do
      guarded = lambda do |x|
        raise "evaluated at endpoint #{x}" if x <= 0.0 || x >= 1.0
        Math.log(x) * Math.log(1.0 - x)
      end
      result = Integration.tanh_sinh(0, 1, &guarded)
      expect(result).to be_converged
      expect(result.value).to be_within(1e-9).of(2.0 - Math::PI**2 / 6.0)
    end

    it 'reports evaluation count, levels and separate error estimates' do
      result = Integration.tanh_sinh(0, 1) { |x| Math.log(x) }
      expect(result.evaluations).to be > 0
      expect(result.levels).to be >= 1
      expect(result.error).to be <= 1e-10
      expect(result.tail_error).to be <= 1e-10
    end

    it 'checks the tail truncation independently of mesh agreement' do
      # A tiny max_points cap forces premature tail truncation: two mesh
      # levels can agree while the neglected tail is still large, so the
      # result must not report convergence.
      result = Integration.tanh_sinh(0, 1, max_points: 4) do |x|
        1.0 / Math.sqrt(x)
      end
      expect(result).not_to be_converged
      expect(result.tail_error).to be > 1e-10
    end

    it 'does not report divergent integrals as tolerance achieved' do
      result = Integration.tanh_sinh(0, 1) { |x| 1.0 / x }
      expect(result).not_to be_converged
      expect(result.status).not_to eq(:ok)
    end

    it 'detects blow-up divergent integrals' do
      result = Integration.tanh_sinh(0, 1) { |x| 1.0 / x**2 }
      expect(result).not_to be_converged
      expect(result.status).to eq(:divergent)
    end

    it 'reports :max_levels when the refinement limit is hit first' do
      result = Integration.tanh_sinh(0, 1, max_levels: 1) do |x|
        Math.exp(-x**2 / 2.0)
      end
      expect(result.status).to eq(:max_levels)
      expect(result).not_to be_converged
    end

    it 'respects reversed bounds' do
      result = Integration.tanh_sinh(1, 0) { |x| Math.log(x) }
      expect(result.value).to be_within(1e-9).of(1.0)
    end

    it 'rejects infinite bounds and zero-length intervals' do
      expect { Integration.tanh_sinh(0, Integration::Infinity) { |x| x } }.to raise_exception
      expect { Integration.tanh_sinh(1, 1) { |x| x } }.to raise_exception
    end
  end

  describe '.integrate with method: :tanh_sinh' do
    it 'returns the value for integrable endpoint singularities' do
      value = Integration.integrate(0, 1, method: :tanh_sinh) { |x| Math.log(x) }
      expect(value).to be_within(1e-9).of(-1.0)
      value = Integration.integrate(0, 1, method: :tanh_sinh) { |x| 1.0 / Math.sqrt(x) }
      expect(value).to be_within(1e-9).of(2.0)
    end

    it 'raises instead of returning an unverified value for divergent integrals' do
      expect do
        Integration.integrate(0, 1, method: :tanh_sinh) { |x| 1.0 / x }
      end.to raise_exception(/tanh_sinh/)
    end

    it 'rejects infinite bounds' do
      expect do
        Integration.integrate(0, Integration::Infinity, method: :tanh_sinh) { |x| Math.exp(-x) }
      end.to raise_exception
    end
  end

  describe '.integrate defaults' do
    it 'keeps the legacy default method unchanged' do
      # Without GSL the default is :simpson; tanh-sinh is strictly opt-in.
      value = Integration.integrate(1, 2, tolerance: 1e-8) { |x| x**2 }
      expect(value).to be_within(1e-6).of(7.0 / 3)
    end
  end
end
