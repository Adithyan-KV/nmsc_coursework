import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def main():
    A = np.array([[3, 2],
                  [2, 6]])
    x_estimate = np.array([[-2], [-2]])
    b = np.array([[2], [-8]])
    x_final, errors, x_values = solve_msd(A, x_estimate, b, 0.001)

    # plot error
    print(f'x value:\n{x_final}')
    plt.plot(errors, 'o-')
    plt.title('Error vs number of iterations')
    plt.xlabel('iteration')
    plt.ylabel('Error')
    plt.show()
    plot_objective_function(A, b, x_values)


def solve_msd(A, x, b, error_threshold):
    max_iterations = 1000

    r = b - np.matmul(A, x)
    error = float('inf')
    errors = []
    x_values = []

    for i in range(max_iterations):
        # if the error is less than threshold return results
        if error < error_threshold:
            print(f'Converged in {i} iterations')
            return(x, errors, x_values)
        q = np.matmul(A, r)
        alpha = np.matmul(np.transpose(r), r) / np.matmul(np.transpose(r), q)
        x_new = x + alpha * r
        error = np.sum((x_new - x)**2)
        x = x_new

        if (i % 50 == 0):
            r = b - np.matmul(A, x)
        else:
            r = r - alpha * q

        # for plotting purposes
        x_values.append(x_new)
        errors.append(error)
    raise Exception("Maximum iterations exceeded with no convergence")


def plot_objective_function(A, b, x_values):
    x = np.arange(-10, 10, 0.1)
    y = np.arange(-10, 10, 0.1)
    X, Y = np.meshgrid(x, y)
    f_values = np.zeros((len(x), len(y)))
    x_descent = np.zeros(len(x_values))
    y_descent = np.zeros_like(x_descent)
    f_descent = np.zeros_like(x_descent)
    for i in range(len(x)):
        for j in range(len(y)):
            x_vec = np.array([x[i], y[j]])
            f = get_objective_function_value(A, b, x_vec)
            # print(f)
            f_values[i, j] = f
    for index, x_value in enumerate(x_values):
        x_descent[index] = x_value[0]
        y_descent[index] = x_value[1]
        f_descent[index] = get_objective_function_value(A, b, x_value)
    fig = plt.figure()
    ax = fig.add_subplot(121, projection='3d')
    ax.set_title('Objective function plot')
    # surface plot
    surface = ax.plot_surface(X, Y, f_values, cmap='coolwarm')
    fig.colorbar(surface, shrink=0.5, aspect=10)
    # plt.show()

    # contour plot with descent path for better visualization
    ax_2 = fig.add_subplot(122, projection='3d')
    ax_2.set_title('Contour plot of objective function and descent path')
    contour = ax_2.contour(X, Y, f_values, cmap='coolwarm')
    ax_2.plot(x_descent, y_descent, f_descent, color='red')
    ax_2.scatter(x_descent, y_descent, f_descent, color='red')
    fig.colorbar(contour, shrink=0.5, aspect=10)
    plt.show()


def get_objective_function_value(A, b, x):
    term_1 = np.matmul(np.transpose(x), np.matmul(A, x))
    term_2 = np.matmul(np.transpose(b), x)
    f = 1 / 2 * term_1 + term_2
    return f


if __name__ == "__main__":
    main()
