$LOAD_PATH.unshift File.expand_path('../spec', __FILE__)

require 'spec_helper'

describe 'Integration.mlmc' do
  # Toy model with known analytic expectations: for omega ~ U(0, 1),
  # the level-l approximation of exp(omega) rounds omega down to a grid
  # of 2^l cells:
  #   P_l(omega) = exp(floor(omega * 2^l) / 2^l)
  #   E[P_l] = (e - 1) / (2^l * (exp(2^-l) - 1))  ->  e - 1
  # Fine and coarse evaluations share the same omega, so the level
  # difference P_l - P_{l-1} has a variance that decays like 4^-l.
  def analytic_level_expectation(level)
    (Math::E - 1) / ((2.0**level) * (Math.exp(2.0**-level) - 1))
  end

  coupled_sampler = lambda do |level, rng|
    omega = rng.rand
    fine = Math.exp((omega * 2**level).floor / (2.0**level))
    if level.zero?
      fine
    else
      coarse = Math.exp((omega * 2**(level - 1)).floor / (2.0**(level - 1)))
      [fine, coarse]
    end
  end

  # Same marginals, but the fine and coarse evaluations use independent
  # random inputs, deliberately breaking the coupling.
  decoupled_sampler = lambda do |level, rng|
    fine = Math.exp((rng.rand * 2**level).floor / (2.0**level))
    if level.zero?
      fine
    else
      coarse = Math.exp((rng.rand * 2**(level - 1)).floor / (2.0**(level - 1)))
      [fine, coarse]
    end
  end

  it 'requires a sampler block' do
    expect { Integration.mlmc(tolerance: 0.1) }.to raise_error(/No sampler/)
  end

  it 'requires exactly one of :tolerance or :budget' do
    expect { Integration.mlmc(&coupled_sampler) }.to raise_error(/tolerance.*budget|budget/)
    expect do
      Integration.mlmc(tolerance: 0.1, budget: 100, &coupled_sampler)
    end.to raise_error(/not both/)
  end

  it 'rejects a sampler that does not return a coupled pair' do
    bad = ->(_level, _rng) { 1.0 }
    expect do
      Integration.mlmc(tolerance: 0.5, min_levels: 2, max_levels: 2, &bad)
    end.to raise_error(/coupled|fine_value/)
  end

  it 'estimates the known analytic expectation within the sampling error' do
    res = Integration.mlmc(tolerance: 0.02, rng: Random.new(42), &coupled_sampler)
    exact = Math::E - 1
    expect(res.estimate).to be_within(0.1).of(exact)
    expect(res.standard_error).to be > 0
    expect(res.standard_error).to be < 0.05
  end

  it 'reports per-level sample counts, means, variances and costs' do
    res = Integration.mlmc(tolerance: 0.02, rng: Random.new(42), &coupled_sampler)
    expect(res.samples.keys).to eq((0...res.levels).to_a)
    res.samples.each_value { |n| expect(n).to be >= 25 }
    expect(res.level_means.keys).to eq(res.samples.keys)
    expect(res.level_variances.keys).to eq(res.samples.keys)
    expect(res.level_costs.keys).to eq(res.samples.keys)
    expect(res.estimate).to be_within(1e-12).of(res.level_means.values.inject(0, :+))
    expect(res.variance).to be_within(1e-12).of(res.standard_error**2)
    expect(res.total_samples).to eq(res.samples.values.inject(0, :+))
  end

  it 'is deterministic for a fixed random source' do
    a = Integration.mlmc(tolerance: 0.02, rng: Random.new(123), &coupled_sampler)
    b = Integration.mlmc(tolerance: 0.02, rng: Random.new(123), &coupled_sampler)
    expect(a.estimate).to eq(b.estimate)
    expect(a.samples).to eq(b.samples)
  end

  it 'shows much larger level-difference variance when coupling is broken' do
    opts = { budget: 1, pilot_samples: 500, min_levels: 3, max_levels: 3 }
    coupled = Integration.mlmc(opts.merge(rng: Random.new(99)), &coupled_sampler)
    decoupled = Integration.mlmc(opts.merge(rng: Random.new(99)), &decoupled_sampler)
    expect(coupled.level_variances[2] * 5).to be < decoupled.level_variances[2]
  end

  it 'allocates more samples to cheap levels under a cost budget' do
    res = Integration.mlmc(
      budget: 20_000, costs: ->(l) { 4.0**l }, rng: Random.new(7), &coupled_sampler
    )
    expect(res.samples[1]).to be > res.samples[2]
    expect(res.samples[2]).to be > res.samples[3]
  end

  it 'extends the level hierarchy when the finest correction has not decayed' do
    coarse = Integration.mlmc(
      tolerance: 0.5, rng: Random.new(5), max_levels: 12, &coupled_sampler
    )
    tight = Integration.mlmc(
      tolerance: 0.005, rng: Random.new(5), max_levels: 12, &coupled_sampler
    )
    expect(tight.levels).to be > coarse.levels
    expect(coarse.levels).to eq(3)
  end

  it 'respects max_levels' do
    res = Integration.mlmc(
      tolerance: 0.01, rng: Random.new(5), max_levels: 4, &coupled_sampler
    )
    expect(res.levels).to eq(4)
  end

  it 'matches the analytic expectation of the finest level it reaches' do
    res = Integration.mlmc(
      tolerance: 0.05, rng: Random.new(11), min_levels: 4, max_levels: 4, &coupled_sampler
    )
    exact = analytic_level_expectation(3)
    expect(res.estimate).to be_within(6 * res.standard_error).of(exact)
  end

  it 'returns a sampling-error confidence interval around the estimate' do
    res = Integration.mlmc(tolerance: 0.02, rng: Random.new(42), &coupled_sampler)
    lo, hi = res.confidence_interval
    expect(lo).to be < res.estimate
    expect(hi).to be > res.estimate
    expect(hi - lo).to be_within(1e-9).of(2 * 1.959963986130195 * res.standard_error)
  end

  it 'exposes a bias indicator separate from the sampling error' do
    res = Integration.mlmc(tolerance: 0.02, rng: Random.new(42), &coupled_sampler)
    expect(res.bias_estimate).to be >= 0
    expect(res.to_h[:bias_estimate]).to eq(res.bias_estimate)
  end
end
