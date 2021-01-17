import math


def main():
    a = 1.0
    b = 3.0
    tolerance = 1e-8

    global function_call_count
    function_call_count = 0

    f = lambda x: (100 / (x**2)) * math.sin(10 / x)

    I = integrate_by_adaptive_quadrature(f, a, b, tolerance)

    print(f'The value of the integral is {I}')
    print(f'Number of function calls for simpsons rule:{function_call_count}')


def integrate_by_adaptive_quadrature(function, lower_lim, upper_lim, tolerance):
    S_tot = integrate_by_simpsons(function, lower_lim, upper_lim)

    I = integrate_two_halves(
        function, lower_lim, upper_lim, S_tot, tolerance)

    return I


def integrate_by_simpsons(function, lower_limit, upper_limit):

    global function_call_count
    function_call_count += 1

    fa = function(lower_limit)
    fb = function(upper_limit)
    midpoint = (lower_limit + upper_limit) / 2
    fm = function(midpoint)
    h = abs(upper_limit - lower_limit) / 2

    S = h * (fa + 4 * fm + fb) / 3

    return S


def integrate_two_halves(function, lower_limit, upper_limit, S_previous, tolerance):

    midpoint = (lower_limit + upper_limit) / 2
    S_left = integrate_by_simpsons(function, lower_limit, midpoint)
    S_right = integrate_by_simpsons(function, midpoint, upper_limit)

    error = S_left + S_right - S_previous

    if abs(error) < 10 * tolerance:
        return S_left + S_right + error / 15
    else:
        return (integrate_two_halves(function, lower_limit, midpoint,
                                     S_left, tolerance / 2)
                + integrate_two_halves(function, midpoint, upper_limit,
                                       S_right, tolerance / 2))


if __name__ == "__main__":
    main()
