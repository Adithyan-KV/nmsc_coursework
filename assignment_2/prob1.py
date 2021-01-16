import numpy as np
import matplotlib.pyplot as plt


def main():
    A = np.array([[3, 2],
                  [2, 6]])
    x_estimate = np.array([[-2], [-2]])
    b = np.array([[2], [-8]])
    x_final, errors = solve_msd(A, x_estimate, b, 0.001)
    print(x_final)
    plt.plot(errors)
    plt.show()


def solve_msd(A, x, b, error_threshold):
    max_iterations = 1000

    r = b - np.matmul(A, x)
    delta = np.matmul(np.transpose(r), r)
    errors = []

    for i in range(max_iterations):
        if delta[0, 0] < error_threshold:
            print(f'Converged in {i} iterations')
            return(x, errors)
        q = np.matmul(A, r)
        alpha = delta / np.matmul(np.transpose(r), q)
        x = x + alpha * r
        if i % 50 == 0:
            r = b - np.matmul(A, x)
        else:
            r = r - alpha * q
        delta = np.matmul(np.transpose(r), r)
        errors.append(delta[0, 0])
    raise Exception("Maximum iterations exceeded with no convergence")


if __name__ == "__main__":
    main()
