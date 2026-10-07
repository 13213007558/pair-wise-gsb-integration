# Integration

[![Build Status](https://travis-ci.org/SciRuby/integration.svg?branch=master)](https://travis-ci.org/SciRuby/integration)
[![Code Climate](https://codeclimate.com/github/SciRuby/integration/badges/gpa.svg)](https://codeclimate.com/github/SciRuby/integration)

Numerical integration for Ruby, with a simple interface. Several different integration methods are implemented.

## Installation

```bash
$ gem install integration
```

Pure Ruby supports infinite bounds when a tail bound is supplied. You can also install GSL (please use your package manager) for a native adaptive implementation:

```bash
$ gem install rb-gsl
```

## Documentation

You can see the latest documentation for Integration at Rubydoc: [rubydoc.info/github/sciruby/integration](http://www.rubydoc.info/github/sciruby/integration).

## Features

* Use only one method: Integration.integrate
* Rectangular, Trapezoidal, Simpson, Adaptive quadrature, Monte Carlo and Romberg integration methods available on pure Ruby.
* Pure-Ruby infinite intervals report the computed interval error, tail bound, and unproven remainder separately.
* If available, uses Ruby/GSL QNG non-adaptive Gauss-Kronrod integration and QAG adaptive integration, with support for Infinity+ and Infinity-.

## Examples

```ruby
# Integrate x^2 in [1, 2] using Simpson's rule.
Integration.integrate(1, 2, tolerance: 1e-10, method: :simpson) do |x|
  x**2
end
# => 2.333333

# Infinite bounds in pure Ruby. Provide an upper bound on the
# absolute value of the remaining tail integral.

normal_pdf = lambda do |x|
  (1 / Math.sqrt(2 * Math::PI)) * Math.exp(-0.5 * (x**2))
end

tail_bound = lambda do |x, direction|
  normal_pdf.call(x) / x.abs
end

result = Integration.integrate_infinite(
  Integration::MInfinity,
  Integration::Infinity,
  tolerance: 1e-10,
  tail_bound: tail_bound,
  &normal_pdf
)

result.value           # => approximately 1.0
result.interval_error  # => error estimate for the computed intervals
result.tail_bound      # => caller-provided bound on the omitted tails
result.unproven        # => 0 when the total bound is below tolerance

# Without a tail bound, only the first finite part is reported and the
# unproven remainder is Infinity; a large fixed endpoint is not substituted.
partial = Integration.integrate_infinite(1, Integration::Infinity) { |x| 1.0 / x }
partial.tail_bound # => Infinity
partial.unproven   # => Infinity

# When rb-gsl is installed, Integration.integrate continues to use QAG
# directly for infinite bounds.
```

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
