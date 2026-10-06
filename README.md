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
* Validated integration via Integration.validated_integrate: rigorous lower/upper bounds that are guaranteed to contain the true integral.

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


## Validated integration

Some tasks need bounds that are *guaranteed* to contain the true value, not
just a point estimate with a sampling error estimate. For those, use
`Integration.validated_integrate`. Instead of a pointwise function you supply
an *enclosure oracle*: a block that receives a subinterval `(a, b)` and
returns a pair `[lo, hi]` such that `lo <= f(x) <= hi` for every `x` in
`[a, b]`. The algorithm adaptively bisects the subinterval with the widest
contribution and accumulates rigorous bounds:

```ruby
# Integral of x^2 over [0, 1] is 1/3. x^2 is increasing on [0, 1],
# so [a^2, b^2] is a valid enclosure on any subinterval [a, b].
result = Integration.validated_integrate(Rational(0), Rational(1),
                                         tolerance: Rational(1, 1000)) do |a, b|
  [a * a, b * b]
end

result.include?(Rational(1, 3)) # => true, guaranteed
result.width                    # => <= 1/1000
result.status                   # => :converged
```

Sampling a plain Float callback cannot prove anything about the unsampled
error, so validated integration never pretends otherwise:

* With Rational/Integer bounds and enclosures (`arithmetic: :exact`, picked
  automatically by `arithmetic: :auto`), all arithmetic is exact.
* With Floats (`arithmetic: :float`), every operation is rounded outward
  through a running error bound, so the enclosure still holds. The oracle
  itself must round its own computations outward, e.g.
  `[Math.exp(a) * (1 - 1e-13), Math.exp(b) * (1 + 1e-13)]` for exp on [a, b].

The evaluation budget is controlled with `max_evaluations` (default 10_000).
If the budget runs out before `tolerance` is met, the current - still valid -
enclosure is returned with `result.status == :budget_exhausted`.

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
