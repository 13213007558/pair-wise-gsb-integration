module Integration
  # Result returned by oscillatory quadrature routines.
  OscillatoryResult = Struct.new(:value, :error_estimate, :amplitude_evaluations) do
    alias_method :evaluations, :amplitude_evaluations
    alias_method :function_evaluations, :amplitude_evaluations

    def to_f
      value
    end

    def to_h
      {
        value: value,
        error_estimate: error_estimate,
        amplitude_evaluations: amplitude_evaluations
      }
    end
  end

  GK15_NODES = [
    -0.9914553711208126, -0.9491079123427585, -0.8648644233597691,
    -0.7415311855993945, -0.5860872354676911, -0.4058451513773972,
    -0.20778495500789848, 0.0, 0.20778495500789848,
    0.4058451513773972, 0.5860872354676911, 0.7415311855993945,
    0.8648644233597691, 0.9491079123427585, 0.9914553711208126
  ]

  GK15_WEIGHTS = [
    0.022935322010529224, 0.06309209262997856, 0.10479001032225019,
    0.14065325971552592, 0.1690047266392679, 0.19035057806478542,
    0.20443294007529889, 0.20948214108472782, 0.20443294007529889,
    0.19035057806478542, 0.1690047266392679, 0.14065325971552592,
    0.10479001032225019, 0.06309209262997856, 0.022935322010529224
  ]

  GK7_WEIGHTS = {
    1 => 0.1294849661688697,
    3 => 0.27970539148927667,
    5 => 0.38183005050511894,
    7 => 0.4179591836734694,
    9 => 0.38183005050511894,
    11 => 0.27970539148927667,
    13 => 0.1294849661688697
  }

  class << self
    # Integrate f(x)*cos(omega*x) or f(x)*sin(omega*x) on a finite interval.
    #
    # Options:
    # [:kind]       +:cosine+ (or +:cos+) or +:sine+ (or +:sin+).
    # [:tolerance]  Requested absolute error. Default: 1e-10.
    # [:panels]     Initial Filon panels. Default: 4.
    # [:max_panels] Maximum panels used during refinement. Default: 4096.
    #
    # The number of panels is increased only to reduce the error estimate;
    # it is never chosen as a multiple of the frequency.
    def oscillatory_integral(lower_bound, upper_bound, frequency, options = {}, &amplitude)
      fail 'No amplitude passed' unless block_given?
      fail 'Non-numeric bounds' unless lower_bound.is_a?(Numeric) && upper_bound.is_a?(Numeric)
      fail 'Non-numeric frequency' unless frequency.is_a?(Numeric)
      fail 'Finite bounds are required' if infinite?(lower_bound) || infinite?(upper_bound)

      orientation = lower_bound <= upper_bound ? 1.0 : -1.0
      a = [lower_bound, upper_bound].min.to_f
      b = [lower_bound, upper_bound].max.to_f
      omega = frequency.to_f
      tolerance = options.fetch(:tolerance, 1e-10).to_f.abs
      panels = options.fetch(:panels, 4).to_i
      max_panels = options.fetch(:max_panels, 4096).to_i
      kind = options.fetch(:kind, :cosine).to_sym

      fail 'Invalid bounds' unless a.finite? && b.finite?
      fail 'Invalid frequency' if !omega.finite?
      fail 'Invalid tolerance' unless tolerance.finite? && tolerance > 0
      fail 'Invalid initial panels' unless panels > 0
      fail 'Invalid maximum panels' unless max_panels >= panels
      fail 'Unknown oscillatory factor' unless [:cosine, :cos, :sine, :sin].include?(kind)

      evaluations = 0
      count_amplitude = proc do |x|
        evaluations += 1
        amplitude.call(x)
      end

      if a == b
        return OscillatoryResult.new(0.0, 0.0, evaluations)
      end

      if (omega * (b - a)).abs <= 1.0
        oscillating_function = case kind
        when :sine, :sin
          proc { |x| count_amplitude.call(x) * Math.sin(omega * x) }
        else
          proc { |x| count_amplitude.call(x) * Math.cos(omega * x) }
        end

        value, error = adaptive_gk15(a, b, tolerance, 0, oscillating_function)
        return OscillatoryResult.new(orientation * value, error, evaluations)
      end

      previous = filon_quadrature(a, b, panels, omega, kind, count_amplitude)

      while panels < max_panels
        panels *= 2
        value = filon_quadrature(a, b, panels, omega, kind, count_amplitude)
        error = (value - previous).abs

        if error <= tolerance
          return OscillatoryResult.new(orientation * value, error, evaluations)
        end

        previous = value
      end

      fail 'Oscillatory integral did not converge'
    end
    alias_method :oscillatory, :oscillatory_integral

    def integrate_oscillatory(lower_bound, upper_bound, frequency, options = {}, &amplitude)
      oscillatory_integral(lower_bound, upper_bound, frequency, options, &amplitude)
    end

    private

    def adaptive_gk15(a, b, tolerance, depth, function)
      value, kronrod_error = gk15_rule(a, b, function)

      return [value, kronrod_error] if kronrod_error <= tolerance || depth >= 100

      middle = (a + b) / 2.0
      left_value, left_error = adaptive_gk15(a, middle, tolerance / 2.0, depth + 1, function)
      right_value, right_error = adaptive_gk15(middle, b, tolerance / 2.0, depth + 1, function)
      [left_value + right_value, left_error + right_error]
    end

    def gk15_rule(a, b, function)
      center = (a + b) / 2.0
      half_width = (b - a) / 2.0
      kronrod = 0.0
      gauss = 0.0

      GK15_NODES.each_with_index do |node, index|
        sample = function.call(center + half_width * node)
        weighted_sample = GK15_WEIGHTS[index] * sample
        kronrod += weighted_sample
        gauss += GK7_WEIGHTS.fetch(index, 0.0) * sample
      end

      value = half_width * kronrod
      difference = (half_width * (kronrod - gauss)).abs
      [value, [difference, value.abs * Float::EPSILON].max]
    end

    def filon_quadrature(a, b, panels, omega, kind, amplitude)
      h = (b - a) / panels.to_f
      values = Array.new(panels + 1) do |index|
        amplitude.call(a + h * index)
      end
      sum = 0.0

      panels.times do |index|
        phase = omega * (a + h * index)
        midpoint_value = amplitude.call(a + h * (index + 0.5))
        cos_weights, sin_weights = filon_quadratic_weights(omega * h)

        left_weight, middle_weight, right_weight = case kind
        when :sine, :sin
          [
            h * (cos_weights[0] * Math.sin(phase) + sin_weights[0] * Math.cos(phase)),
            h * (cos_weights[1] * Math.sin(phase) + sin_weights[1] * Math.cos(phase)),
            h * (cos_weights[2] * Math.sin(phase) + sin_weights[2] * Math.cos(phase))
          ]
        else
          [
            h * (cos_weights[0] * Math.cos(phase) - sin_weights[0] * Math.sin(phase)),
            h * (cos_weights[1] * Math.cos(phase) - sin_weights[1] * Math.sin(phase)),
            h * (cos_weights[2] * Math.cos(phase) - sin_weights[2] * Math.sin(phase))
          ]
        end

        sum += left_weight * values[index] +
               middle_weight * midpoint_value +
               right_weight * values[index + 1]
      end

      sum
    end

    def filon_quadratic_weights(k)
      if k.abs < 1.0
        moments = [
          [
            oscillatory_moment(:cosine, 0, k),
            oscillatory_moment(:cosine, 1, k),
            oscillatory_moment(:cosine, 2, k)
          ],
          [
            oscillatory_moment(:sine, 0, k),
            oscillatory_moment(:sine, 1, k),
            oscillatory_moment(:sine, 2, k)
          ]
        ]
      else
        c0 = Math.sin(k) / k
        s0 = (1.0 - Math.cos(k)) / k
        c1 = (Math.cos(k) + k * Math.sin(k) - 1.0) / (k * k)
        s1 = (Math.sin(k) - k * Math.cos(k)) / (k * k)
        c2 = Math.sin(k) / k + 2.0 * Math.cos(k) / (k * k) -
             2.0 * Math.sin(k) / (k * k * k)
        s2 = -Math.cos(k) / k + 2.0 * Math.sin(k) / (k * k) +
             2.0 * Math.cos(k) / (k * k * k) - 2.0 / (k * k * k)
        moments = [[c0, c1, c2], [s0, s1, s2]]
      end

      c0, c1, c2 = moments[0]
      s0, s1, s2 = moments[1]
      [
        [
          c0 - 3.0 * c1 + 2.0 * c2,
          4.0 * c1 - 4.0 * c2,
          2.0 * c2 - c1
        ],
        [
          s0 - 3.0 * s1 + 2.0 * s2,
          4.0 * s1 - 4.0 * s2,
          2.0 * s2 - s1
        ]
      ]
    end

    def oscillatory_moment(kind, degree, k)
      sum = 0.0
      term_degree = kind == :sine ? 1 : 0
      sign_index = 0
      term = kind == :sine ? k.to_f : 1.0

      8.times do
        sum += (sign_index.even? ? 1.0 : -1.0) * term / (term_degree + degree + 1).to_f
        term *= k * k / ((term_degree + 1).to_f * (term_degree + 2).to_f)
        term_degree += 2
        sign_index += 1
      end

      sum
    end
  end
end
