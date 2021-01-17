import numpy as np
import matplotlib.pyplot as plt


def main():

    # plot velocities vs time
    velocities_scenario_A, t_scenario_A = runge_kutta_4(0, 1.0, 0.7, 0.2, 1e-8)
    velocities_scenario_B, t_scenario_B = runge_kutta_4(
        0, 0.7, 0.18, 0.2, 1e-8)

    fig, plots = plt.subplots(1, 2)
    fig.suptitle('Velocity vs Time')
    plots[0].plot(t_scenario_A, velocities_scenario_A)
    plots[0].set_title('Scenario A')
    plots[0].set(xlabel="Time(s)", ylabel="Velocity(m/s)")
    plots[1].plot(t_scenario_B, velocities_scenario_B)
    plots[1].set_title('Scenario B')
    plots[1].set(xlabel="Time(s)", ylabel="Velocity(m/s)")
    plt.show()

    # get velocities using analytical calculation
    v_analytic_A = compute_velociy_analytic(t_scenario_A, 1.0, 0.7)
    v_analytic_B = compute_velociy_analytic(t_scenario_B, 0.7, 0.18)

    # calculate errors w.r.t analytic
    errors_A = np.abs(v_analytic_A - velocities_scenario_A)
    errors_B = np.abs(v_analytic_B - velocities_scenario_B)

    # plot errors vs time
    fig, plots = plt.subplots(1, 2)
    fig.suptitle('Error w.r.t analytical solution')
    plots[0].plot(t_scenario_A, errors_A)
    plots[0].set_title('Scenario A')
    plots[0].set(xlabel="Time(s)", ylabel="Absolute error")
    plots[1].plot(t_scenario_B, errors_B)
    plots[1].set_title('Scenario B')
    plots[1].set(xlabel="Time(s)", ylabel="Absolute error")
    plt.show()

    v_A_1, t_A_1 = runge_kutta_4(0, 1.0, 0.7, 0.1, 1e-6)
    v_A_2, t_A_2 = runge_kutta_4(0, 1.0, 0.7, 0.2, 1e-6)
    v_A_5, t_A_5 = runge_kutta_4(0, 1.0, 0.7, 0.5, 1e-6)
    v_A_10, t_A_10 = runge_kutta_4(0, 1.0, 0.7, 1.0, 1e-6)

    v_B_1, t_B_1 = runge_kutta_4(0, 0.7, 0.18, 0.1, 1e-6)
    v_B_2, t_B_2 = runge_kutta_4(0, 0.7, 0.18, 0.2, 1e-6)
    v_B_5, t_B_5 = runge_kutta_4(0, 0.7, 0.18, 0.5, 1e-6)
    v_B_10, t_B_10 = runge_kutta_4(0, 0.7, 0.18, 0.1, 1e-6)

    v_analytic_A_1 = compute_velociy_analytic(t_A_1, 1.0, 0.7)
    v_analytic_A_2 = compute_velociy_analytic(t_A_2, 1.0, 0.7)
    v_analytic_A_5 = compute_velociy_analytic(t_A_5, 1.0, 0.7)
    v_analytic_A_10 = compute_velociy_analytic(t_A_10, 1.0, 0.7)

    v_analytic_B_1 = compute_velociy_analytic(t_B_1, 0.7, 0.18)
    v_analytic_B_2 = compute_velociy_analytic(t_B_2, 0.7, 0.18)
    v_analytic_B_5 = compute_velociy_analytic(t_B_5, 0.7, 0.18)
    v_analytic_B_10 = compute_velociy_analytic(t_B_10, 0.7, 0.18)

    RMSE_A = RMSE(v_analytic_A_1, v_A_1), RMSE(v_analytic_A_2, v_A_2), RMSE(
        v_analytic_A_5, v_A_5), RMSE(v_analytic_A_10, v_A_10)
    RMSE_B = RMSE(v_analytic_B_1, v_B_1), RMSE(v_analytic_B_2, v_B_2), RMSE(
        v_analytic_B_5, v_B_5), RMSE(v_analytic_B_10, v_B_10)

    t_steps = [0.1, 0.2, 0.5, 1.0]

    fig, plots = plt.subplots(1, 2)
    fig.suptitle('RMSE vs time')
    plots[0].plot(t_steps, RMSE_A, '-o')
    plots[0].set_title('Scenario A')
    plots[0].set(xlabel="Time(s)", ylabel="RMSE")
    plots[1].plot(t_steps, RMSE_B, '-o')
    plots[1].set_title('Scenario B')
    plots[1].set(xlabel="Time(s)", ylabel="RMSE")
    plt.show()


def f(velocity, viscosity_coefficient, area):
    # All values in SI units
    g = 9.81
    density = 1.21
    mass = 85

    f = g - (viscosity_coefficient * density * area * velocity**2) / (2 * mass)

    return f


def runge_kutta_4(initial_v, visc_coeff, area, h, tolerance):
    max_iterations = 1000
    error = float('inf')
    t = 0
    v = initial_v
    v_list = [v]
    t_list = [t]

    for _ in range(max_iterations):
        if error < tolerance:
            return v_list, t_list
        k1 = h * f(v, visc_coeff, area)
        k2 = h * f(v + k1 / 2, visc_coeff, area)
        k3 = h * f(v + k2 / 2, visc_coeff, area)
        k4 = h * f(v + k3, visc_coeff, area)

        v = v + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t = t + h

        v_list = np.append(v_list, v)
        t_list = np.append(t_list, t)

        error = (abs(v_list[-1] - v_list[-2])) / v_list[-1]
    raise Exception("Failed to reach solution after maximum iterations")


def compute_velociy_analytic(t, visc_coeff, area):
    g = 9.81
    density = 1.21
    mass = 85

    v_t = np.sqrt((2 * mass * g) / (visc_coeff * density * area))
    tau = np.sqrt((2 * mass) / (visc_coeff * density * area * g))

    v = v_t * ((np.exp(2 * t / tau) - 1) / (np.exp(2 * t / tau) + 1))

    return v


def RMSE(array_1, array_2):

    RMSE = (np.sum((array_1 - array_2)**2)) / array_2.size

    return RMSE


if __name__ == "__main__":
    main()
