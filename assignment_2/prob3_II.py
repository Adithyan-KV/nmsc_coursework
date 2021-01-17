import numpy as np
import matplotlib.pyplot as plt


def main():

    dist_arr_A, t_arr_A, time_req_A, _, dist_TOL_A, time_diff_A, v_max_A = RK4_and_integration(
        0, 1.0, 0.7, 0.2, 1e-8)
    dist_arr_B, t_arr_B, time_req_B, _, _, _, v_max_B = RK4_and_integration(
        0, 0.7, 0.18, 0.2, 1e-8)

    print("Scenario A")
    print("-" * 100)
    print(f"Time till 5000 feet                 : {time_req_A} seconds")
    print(f"Distance till                       : {dist_TOL_A/0.3048} ft ")
    print(f"Time travelled at terminal velocity : {time_diff_A} seconds")
    print(f"Maximum velocity attained           : {v_max_A*3.6} km/hr")

    print("Scenario B")
    print("-" * 100)
    print(f"Time till 5000 feet                 : {time_req_B} seconds")
    print(f"Distance till terminal velocity     : Never attains terminal veocity ")
    print(f"Time travelled at terminal velocity : Never attains terminal velocity")
    print(f"Maximum velocity attained           : {v_max_B*3.6} km/hr")

    distance_analytic_A = compute_distance_analytic(t_arr_A, 1.0, 0.7)
    distance_analytic_B = compute_distance_analytic(t_arr_B, 0.7, 0.18)

    errors_A = np.abs(distance_analytic_A - dist_arr_A)
    errors_B = np.abs(distance_analytic_B - dist_arr_B)

    fig, plots = plt.subplots(1, 2)
    fig.suptitle('Error w.r.t analytic')
    plots[0].plot(t_arr_A, errors_A)
    plots[0].set_title('Scenario A')
    plots[0].set(xlabel="Time(s)", ylabel="Error")
    plots[1].plot(t_arr_B, errors_B)
    plots[1].set_title('Scenario B')
    plots[1].set(xlabel="Time(s)", ylabel="Error")
    plt.show()


def integrate_by_gaussian_quadrature(a, b, fa, fb):

    # coefficient values for gaussian quadrature
    x = np.array([-0.774597, 0, 0.774597])
    w = np.array([0.555556, 0.888889, 0.555556])

    # shifted version of x
    x_bar = ((b - a) * x / 2) + (a + b) / 2

    # linear interpolation
    y = (fa * (b - x_bar) + fb * (x_bar - a)) / (b - a)

    # computing the integral
    I = ((b - a) / 2) * np.sum(w * y)

    return I


def f(velocity, visc_coeff, area):
    g = 9.81
    density = 1.21
    mass = 85

    # equation for viscous force
    f = g - (visc_coeff * density * area * velocity**2) / (2 * mass)

    return f


def RK4_and_integration(v0, C, A, h, tolerance):

    max_iterations = 1000
    error = float('inf')
    dist_covered = 0.0
    dist_thresh = (13500 - 5000) * 0.3048
    t = 0
    v = v0
    v_array = np.array([v])
    t_array = np.array([t])
    dist_array = np.array([0])

    dist_lock = True
    TOL_lock = True

    for _ in range(max_iterations):
        if (error > tolerance) or (dist_covered <= dist_thresh):
            k1 = h * f(v, C, A)
            k2 = h * f(v + k1 / 2, C, A)
            k3 = h * f(v + k2 / 2, C, A)
            k4 = h * f(v + k3, C, A)

            v = v + (k1 + 2 * k2 + 2 * k3 + k4) / 6
            t = t + h

            v_array = np.append(v_array, v)
            t_array = np.append(t_array, t)

            error = (abs(v_array[-1] - v_array[-2])) / v_array[-1]
            dist_step = integrate_by_gaussian_quadrature(
                t_array[-2], t_array[-1], v_array[-2], v_array[-1])
            dist_covered += dist_step
            dist_array = np.append(dist_array, dist_covered)

            if (dist_covered > dist_thresh) and dist_lock:
                time_required = t_array[-1]
                dist_lock = False

            if (error <= tolerance) and TOL_lock:
                dist_tolerance = dist_covered
                time_tolerance = t_array[-1]
                TOL_lock = False
        else:
            break

    time_diff = time_required - time_tolerance
    if time_diff > 0:
        v_max = v_array[-1]
    else:
        v_max = v_array[np.nonzero(t_array == time_required)[0][0]]

    return dist_array, t_array, time_required, time_tolerance, dist_tolerance, time_diff, v_max


def compute_distance_analytic(t, visc_coeff, area):
    g = 9.81
    density = 1.21
    mass = 85

    v_t = np.sqrt((2 * mass * g) / (visc_coeff * density * area))
    tau = np.sqrt((2 * mass) / (visc_coeff * density * area * g))
    distance = tau * v_t * np.log(np.cosh(t / tau))

    return distance


if __name__ == "__main__":
    main()
