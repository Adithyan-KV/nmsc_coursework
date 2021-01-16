import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def main():
    A = np.array([[3, 2],
                  [2, 6]])
    x_estimate = np.array([[-2], [-2]])
    b = np.array([[2], [-8]])
    x_final, errors = solve_msd(A, x_estimate, b, 0.001)
    print(x_final)
    plt.plot(errors)
    plt.show()
    plot_objective_function(A, b)


def solve_msd(A, x, b, error_threshold):
    max_iterations = 1000

    r = b - np.matmul(A, x)
    delta = np.matmul(np.transpose(r), r)
    errors = []

    for i in range(max_iterations):
        # if the error is less than threshold return results
        if delta[0, 0] < error_threshold:
            print(f'Converged in {i} iterations')
            return(x, errors)
        q = np.matmul(A, r)
        alpha = delta / np.matmul(np.transpose(r), q)
        x = x + alpha * r

        # recompute the residual every 50 iterations to avoid numerical error
        if i % 50 == 0:
            r = b - np.matmul(A, x)
        else:
            r = r - alpha * q
        delta = np.matmul(np.transpose(r), r)
        errors.append(delta[0, 0])
    raise Exception("Maximum iterations exceeded with no convergence")


def plot_objective_function(A, b):
    x = np.arange(-10, 10, 0.1)
    y = np.arange(-10, 10, 0.1)
    X, Y = np.meshgrid(x, y)
    f_values = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            x_vec = np.array([x[i], y[j]])
            f = get_objective_function_value(A, b, x_vec)
            # print(f)
            f_values[i, j] = f
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surface = ax.plot_surface(X, Y, f_values, cmap='rainbow')
    fig.colorbar(surface, shrink=0.5, aspect=10)
    plt.show()


def get_objective_function_value(A, b, x):
    term_1 = np.matmul(np.transpose(x), np.matmul(A, x))
    term_2 = np.matmul(np.transpose(x), x)
    f = 1 / 2 * term_1 + term_2
    return f


if __name__ == "__main__":
    main()
