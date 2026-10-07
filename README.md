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

When the same simulation exists in a cheap coarse-grid version and an
expensive fine-grid version, `Integration.mlmc` estimates the expectation
of the finest model with a multilevel Monte Carlo strategy. The caller
provides a *coupled* sampler for the level differences
`Y_l = P_l(omega) - P_{l-1}(omega)`, where coarse and fine evaluations
share the same random input `omega`:

```ruby
result = Integration.mlmc(
  sampler: ->(level, n, rng) {
    Array.new(n) do
      omega = rng.rand
      fine  = simulate(level, omega)          # expensive on fine levels
      level.zero? ? fine : fine - simulate(level - 1, omega)
    end
  },
  costs:  ->(level) { 4.0**level },  # cost per coupled sample
  budget: 100_000,                   # total cost units
  pilot_samples: 25,                 # pilot batch per level
  tolerance: 1e-3,                   # optional target RMSE
  rng: Random.new(42)                # fixed seed for reproducibility
)

result.estimate         # MLMC estimate of the expectation
result.standard_error   # Monte Carlo sampling error only (see below)
result.samples          # {0 => 5200, 1 => 1300, ...} samples spent per level
result.confidence_interval(0.95)
```

A small pilot batch is drawn on every level first; the remaining budget is
then distributed with the variance-optimal allocation
`N_l ~ sqrt(Var(Y_l) / C_l)`. Finer levels are appended while the bias
remaining after the finest level, extrapolated from the geometric decay of
the level-difference means, exceeds the allowed bias (`tolerance/sqrt(2)`,
or the current sampling error when no tolerance is given).

Note: the reported standard error and confidence interval quantify Monte
Carlo *sampling* error only. They are not a strict guarantee covering the
discretization bias of the finest level; `result.bias_estimate` gives a
heuristic estimate of that bias.

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
