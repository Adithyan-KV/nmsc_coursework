def main():
    f_1 = lambda p: p * (1 + (7 - p**5) / (p**2))**3
    f_2 = lambda p: p - ((p**5 - 7) / (p**2))**3
    f_3 = lambda p: p - ((p**5 - 7) / (5 * p**4))
    f_4 = lambda p: p - ((p**5 - 7) / (12))
    functions = [f_1, f_2, f_3, f_4]
    for function in functions:
        try:
            iterate_till_convergence(function)
        except:
            print("fucking shit not working")


def iterate_till_convergence(function, estimate=1, tolerance=1e-8, max_iter=10000):
    p = estimate
    error = float('inf')
    for i in range(max_iter):
        if error < tolerance:
            print(f'converged to {p_new} after {i} iterations')
            return i
        p_new = function(p)
        error = abs(p_new - p)
        p = p_new
    raise Exception(
        f"Couldn't reach convergence even after {max_iter} iterations")


if __name__ == "__main__":
    main()
