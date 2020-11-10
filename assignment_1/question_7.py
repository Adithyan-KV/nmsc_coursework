import math


def main():
    # the functions and their derivatives
    f_1 = lambda x: x + math.exp(-x**2) * math.cos(x)
    df_1 = lambda x: 1 - math.exp(-x**2) * \
        (math.sin(x) + 2 * x * math.cos(x))
    ddf_1 = lambda x: math.exp(-x**2) * ((4 * x**2 - 3)
                                         * math.cos(x) + 4 * x * math.sin(x))

    f_2 = lambda x: f_1(x)**2
    df_2 = lambda x: 2 * f_1(x) * df_1(x)
    ddf_2 = lambda x: 2 * (df_1(x)**2 + f_1(x) * ddf_1(x))

    # calculating root by newtons method
    newtons_method(f_1, df_1)
    newtons_method(f_2, df_2)

    # calculating root by secant method
    secant_method(f_1)
    secant_method(f_2)

    # calculating root be modified newtons method
    modified_newtons_method(f_1, df_1, ddf_1)
    modified_newtons_method(f_2, df_2, ddf_2)


def newtons_method(function, derivative, estimate=0, max_iter=10000, tolerance=1e-8):
    p_0 = estimate
    error = float('inf')
    for i in range(max_iter):
        p = p_0 - function(p_0) / derivative(p_0)
        error = abs(p - p_0)
        if error < tolerance:
            print(f"converged to {p_0} after {i} iterations")
            return(p_0)
        p_0 = p
    raise Exception(
        f"Couldn't reach convergence after running for {max_iter} iterations")


def secant_method(function, p_0=0, p_1=1, max_iter=10000, tolerance=1e-8):
    q_0 = function(p_0)
    q_1 = function(p_1)
    for i in range(2, max_iter):
        p = p_1 - q_1 * ((p_1 - p_0) / (q_1 - q_0))
        error = abs(p - p_1)
        if error < tolerance:
            print(f"converged to {p} after {i} iterations")
            return(p)
        p_0 = p_1
        q_0 = q_1
        p_1 = p
        q_1 = function(p)
    raise Exception(
        f"Couldn't reach convergence after running for {max_iter} iterations")


def modified_newtons_method(function, derivative, second_der, p_0=1, max_iter=10000, tolerance=1e-8):
    error = float('inf')
    for i in range(max_iter):
        p = p_0 - (function(p_0) * derivative(p_0)) / \
            (derivative(p_0)**2 - function(p_0) * second_der(p_0))
        error = abs(p - p_0)
        if error < tolerance:
            print(f"converged to {p_0} after {i} iterations")
            return(p_0)
        p_0 = p
    raise Exception(
        f"Couldn't reach convergence after running for {max_iter} iterations")


if __name__ == "__main__":
    main()
