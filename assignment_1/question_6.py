def main():
    # defining the functions for iteratively calculating terms of p
    f_1 = lambda p: p * (1 + (7 - p**5) / (p**2))**3
    f_2 = lambda p: p - ((p**5 - 7) / (p**2))
    f_3 = lambda p: p - ((p**5 - 7) / (5 * p**4))
    f_4 = lambda p: p - ((p**5 - 7) / (12))

    # run all functions till convergence
    functions = [f_1, f_2, f_3, f_4]
    for function in functions:
        iterate_till_convergence(function)


def iterate_till_convergence(function, estimate=1, tolerance=1e-8, max_iter=10000):
    p = estimate
    error = float('inf')
    for i in range(max_iter):
        # if the series has converged return the value and number of iterations
        if error < tolerance:
            print(f'converged to {p_new} after {i} iterations')
            return (p, i)

        # if the value goes out of numeric range, series not converging
        try:
            p_new = function(p)
        except OverflowError:
            print(
                f"Couldn't reach convergence. Value out of numeric range at p{i}")
            return None

        # metric for checking convergence
        error = abs(p_new - p)
        p = p_new

    # if convergence couldn't be reached after specified maximum iterations
    raise Exception(
        f"Couldn't reach convergence even after {max_iter} iterations")


if __name__ == "__main__":
    main()
