require 'spec_helper'

# Toy model with a known analytic expectation, mimicking a simulation with
# cheap coarse-grid and expensive fine-grid versions:
#
#   P_l(omega) = omega * (1 + 2^-l),  omega ~ Uniform(0, 1)
#
# The exact expectation of the continuous model is E[omega] = 0.5, while the
# level-L discretization has E[P_L] = (1 + 2^-L) / 2, so the remaining
# discretization bias after level L is exactly 2^-(L+1).
describe Integration do
  describe '.mlmc' do
    let(:pilot) { 20 }
    let(:costs) { ->(level) { 4.0**level } }

    # Coupled sampler: fine and coarse approximations share the same omega.
    let(:coupled_sampler) do
      lambda do |level, n, rng|
        Array.new(n) do
          omega = rng.rand
          fine = omega * (1 + 2.0**-level)
          level.zero? ? fine : fine - omega * (1 + 2.0**-(level - 1))
        end
      end
    end

    # Deliberately decoupled sampler: the coarse evaluation uses an
    # independent random draw, breaking the shared random input.
    let(:decoupled_sampler) do
      lambda do |level, n, rng|
        Array.new(n) do
          fine = rng.rand * (1 + 2.0**-level)
          level.zero? ? fine : fine - rng.rand * (1 + 2.0**-(level - 1))
        end
      end
    end

    def run_mlmc(sampler, overrides = {})
      defaults = {
        sampler: sampler,
        costs: costs,
        budget: 100_000,
        pilot_samples: pilot,
        levels: [0, 1, 2, 3, 4],
        extend_levels: false,
        rng: Random.new(1234)
      }
      Integration.mlmc(defaults.merge(overrides))
    end

    it 'estimates the analytic expectation of the finest level' do
      result = run_mlmc(coupled_sampler)
      analytic = 0.5 * (1 + 2.0**-4)
      expect((result.estimate - analytic).abs).to be < 4 * result.standard_error
    end

    it 'is reproducible for a fixed random source' do
      first = run_mlmc(coupled_sampler)
      second = run_mlmc(coupled_sampler)
      expect(first.estimate).to eq second.estimate
      expect(first.samples).to eq second.samples
    end

    it 'reports the number of samples spent on each level' do
      result = run_mlmc(coupled_sampler)
      expect(result.levels).to eq [0, 1, 2, 3, 4]
      result.levels.each do |level|
        expect(result.samples[level]).to be >= pilot
      end
      # Cheap coarse levels should receive most of the budget.
      expect(result.samples[0]).to be > result.samples[4]
    end

    it 'stays within the cost budget' do
      result = run_mlmc(coupled_sampler)
      spent = result.samples.inject(0.0) { |acc, (l, n)| acc + n * costs.call(l) }
      expect(spent).to be <= 100_000
      expect(result.cost_spent).to be_within(1e-9).of(spent)
    end

    it 'allocates more samples where variance per cost is high' do
      result = run_mlmc(coupled_sampler)
      # Var(Y_l) = 2^-2l / 12 and C_l = 4^l, so sqrt(V_l / C_l) = 2^-2l / sqrt(12):
      # the allocation should decay roughly by a factor 4 per level.
      ratio = result.samples[0].to_f / result.samples[2]
      expect(ratio).to be > 4
    end

    context 'with a deliberately decoupled sampler' do
      it 'exhibits much larger level-difference variances' do
        coupled = run_mlmc(coupled_sampler)
        decoupled = run_mlmc(decoupled_sampler)
        [1, 2, 3].each do |level|
          expect(decoupled.level_variances[level]).to be > 4 * coupled.level_variances[level]
        end
        expect(decoupled.standard_error).to be > coupled.standard_error
      end
    end

    context 'level extension' do
      it 'adds finer levels while the bias proxy exceeds the allowed bias' do
        result = Integration.mlmc(
          sampler: coupled_sampler,
          costs: ->(level) { 2.0**level },
          budget: 1_000_000,
          pilot_samples: 20,
          levels: [0, 1],
          tolerance: 1e-3,
          max_level: 12,
          rng: Random.new(99)
        )
        # bias proxy after level L is ~2^-(L+1); tolerance/sqrt(2) ~ 7.1e-4
        # requires 2^-(L+1) <= 7.1e-4, i.e. L >= 10.
        expect(result.levels.max).to be > 4
        expect(result.bias_estimate).to be <= 1e-3 / Math.sqrt(2.0)
      end

      it 'does not add levels when extend_levels is false' do
        result = run_mlmc(coupled_sampler, levels: [0, 1])
        expect(result.levels).to eq [0, 1]
      end

      it 'respects max_level' do
        result = Integration.mlmc(
          sampler: coupled_sampler,
          costs: costs,
          budget: 1_000_000,
          pilot_samples: 20,
          levels: [0, 1],
          tolerance: 1e-12,
          max_level: 3,
          rng: Random.new(7)
        )
        expect(result.levels.max).to eq 3
      end
    end

    describe 'confidence interval' do
      it 'is centered on the estimate with normal-quantile half width' do
        result = run_mlmc(coupled_sampler)
        lo, hi = result.confidence_interval(0.95)
        expect(lo).to be < result.estimate
        expect(hi).to be > result.estimate
        expect(hi - lo).to be_within(1e-6).of(2 * 1.959963984540054 * result.standard_error)
      end

      it 'covers sampling error only, not discretization bias' do
        # With a coarse finest level and a large budget, the sampling error
        # becomes tiny while the discretization bias E[P_2] - 0.5 = 0.125
        # remains: the interval must NOT be read as containing the exact
        # expectation 0.5 of the continuous model.
        result = run_mlmc(coupled_sampler, levels: [0, 1, 2], budget: 2_000_000)
        lo, hi = result.confidence_interval(0.95)
        expect(lo).to be > 0.5
        expect(result.bias_estimate).to be_within(0.02).of(0.125)
        expect(lo).to be < 0.625
        expect(hi).to be > 0.625
      end
    end

    it 'requires a sampler, costs and a budget' do
      expect { Integration.mlmc(costs: costs, budget: 10) }.to raise_error(/sampler/)
      expect { Integration.mlmc(sampler: coupled_sampler, budget: 10) }.to raise_error(/costs/)
      expect { Integration.mlmc(sampler: coupled_sampler, costs: costs) }.to raise_error(/budget/)
    end
  end
end
