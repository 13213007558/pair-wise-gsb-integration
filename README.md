# Integration

[![Build Status](https://travis-ci.org/SciRuby/integration.svg?branch=master)](https://travis-ci.org/SciRuby/integration)
[![Code Climate](https://codeclimate.com/github/SciRuby/integration/badges/gpa.svg)](https://codeclimate.com/github/SciRuby/integration)

Numerical integration for Ruby, with a simple interface. Several different integration methods are implemented.

## Installation

```bash
$ gem install integration
```

You can install GSL (please use your package manager) for better performance and support for infinite bounds. After successfully installing the library:

```bash
$ gem install rb-gsl
```

## Documentation

You can see the latest documentation for Integration at Rubydoc: [rubydoc.info/github/sciruby/integration](http://www.rubydoc.info/github/sciruby/integration).

## Features

* Use only one method: Integration.integrate
* Rectangular, Trapezoidal, Simpson, Adaptive quadrature, Monte Carlo and Romberg integration methods available on pure Ruby.
* Multilevel Monte Carlo (MLMC) expectation estimation for hierarchies of discretizations (cheap coarse / expensive fine grids) with a coupled sampler, pilot-based optimal sample allocation and per-level sample accounting.
* If available, uses Ruby/GSL QNG non-adaptive Gauss-Kronrod integration and QAG adaptive integration, with support for Infinity+ and Infinity-

## Examples

```ruby
# Integrate x^2 in [1, 2] using Simpson's rule.
Integration.integrate(1, 2, tolerance: 1e-10, method: :simpson) do |x|
  x**2
end
# => 2.333333

# Support for infinity bounds with GSL QAG adaptative integration.

normal_pdf = lambda do |x|
  (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-0.5 * (x**2))
end

Integration.integrate(Integration::MInfinity, 0, tolerance: 1e-10, &normal_pdf) # => 0.5
Integration.integrate(0, Integration::Infinity , tolerance: 1e-10, &normal_pdf) # => 0.5
```

## Multilevel Monte Carlo

When the same simulation is available at several resolutions (a cheap coarse
grid and increasingly expensive fine grids), `Integration.mlmc` estimates the
expectation of the finest discretization from a *coupled* sampler. The sampler
receives a level and a random source, and must evaluate the fine and coarse
approximations at the *same* random input:

```ruby
sampler = lambda do |level, rng|
  omega = rng.rand # shared random input for both resolutions
  fine  = simulate(omega, grid: 2**level)
  if level.zero?
    fine
  else
    coarse = simulate(omega, grid: 2**(level - 1))
    [fine, coarse] # coupled pair; the difference has small variance
  end
end

result = Integration.mlmc(tolerance: 0.01, rng: Random.new(42), &sampler)
result.estimate        # estimated expectation
result.standard_error  # Monte Carlo sampling error only
result.samples         # => {0 => ..., 1 => ..., ...} samples drawn per level
```

The method first draws a small pilot at every level, then allocates the
remaining samples as N_l proportional to sqrt(V_l / C_l) from the pilot
variances V_l and the per-sample costs C_l (pass `:budget` instead of
`:tolerance` to fix the total cost instead of the sampling error, and
`:costs` to describe the cost per level). Levels are appended automatically
while the finest-level correction has not decayed enough — specifically while
the remaining bias extrapolated from the two finest level means exceeds the
target sampling error — up to `:max_levels`.

Note that `standard_error` and `confidence_interval` quantify the sampling
error only; they are not a guarantee that the exact expectation lies inside
the interval, because the discretization bias of the finest level is not
included. Inspect `result.bias_estimate` for a heuristic indicator of that
bias.

## License

Copyright (c) 2005  Beng
              2011-2015  clbustos

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
THE DEVELOPERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF
OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Except as contained in this notice, the name of the Beng shall not
be used in advertising or otherwise to promote the sale, use or other dealings
in this Software without prior written authorization from Beng.
