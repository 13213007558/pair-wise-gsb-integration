# Example: integrating a vector-valued function with one expensive callback.
#
# Run with: ruby -Ilib examples/vector_integration.rb
#
# A single callback computes 12 quantities that share an expensive common
# part. Integrating each component separately repeats that work; the
# vector-valued entry point shares sample points and interval subdivision
# across all components while checking each one against its own error
# target.

require 'integration'

LOWER = 0.0
UPPER = Math::PI
K_MAX = 12

# One expensive callback returning K_MAX quantities at once.
def expensive_callback(x)
  common = (1..50).inject(0.0) { |acc, j| Math.sin(x * j) * Math.cos(x / (j + 1)) }
  (1..K_MAX).map { |k| Math.exp(-k * x / 5.0) + 1e-12 * common }
end

# Reference: integrate each component on its own, counting evaluations.
scalar_evaluations = 0
reference = (1..K_MAX).map do |k|
  Integration.integrate(LOWER, UPPER, tolerance: 1e-9, method: :adaptive_quadrature) do |x|
    scalar_evaluations += 1
    expensive_callback(x)[k - 1]
  end
end

# Vector-valued: one pass, shared sampling, per-component diagnostics.
result = Integration.integrate_vector(LOWER, UPPER, tolerance: 1e-9) do |x|
  expensive_callback(x)
end

puts format('%-4s %-22s %-22s %-12s %s', 'k', 'vector', 'scalar', 'error est.', 'converged')
result.diagnostics.each do |d|
  puts format('%-4d %-22.15e %-22.15e %-12.3e %s',
              d[:index] + 1, d[:value], reference[d[:index]], d[:error], d[:converged])
end

max_diff = result.values.each_with_index.map { |v, i| (v - reference[i]).abs }.max
puts
puts "max |vector - scalar| : #{format('%.3e', max_diff)}"
puts "callback evaluations  : vector=#{result.evaluations} scalar=#{scalar_evaluations} " +     "(#{format('%.1f', scalar_evaluations.to_f / result.evaluations)}x fewer)"
