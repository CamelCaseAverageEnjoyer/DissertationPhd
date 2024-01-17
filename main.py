from main_objects import *
from plot_func import *


if __name__ == '__main__':
    # Инициализация
    start_navigation_tolerance = 0.95
    t_integrate = 1e4
    dt = 10.

    start_navigation = ['perfect', 'near', 'random'][1]
    model_c = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
    method_navigation = ['kalman_filter rv', 'kalman_filter rvq', 'kalman_filter rv all', 'kalman_filter rvq all'][3]
    kalman_coef = {'q': 1e-16, 'p': [1e-8] * 3, 'r': 1e0}

    o = Objects(n_c=1, n_f=5, kalman_coef=kalman_coef,
                model_c=model_c, dt=dt, start_navigation=start_navigation,
                start_navigation_tolerance=start_navigation_tolerance, method_navigation=method_navigation)
    o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas'][0]
    o.c.gain_mode = o.f.gain_mode
    o.p.k.gain_mode = o.f.gain_mode
    o.p.is_aero = False
    o.p.show_rate = 1

    # Интегрирование
    o.p.integrate(t=t_integrate)

    # Результаты
    print(f"Математическое ожидание: {np.array(o.f.line_difference[0]).mean()}, "
          f"дисперсия: {np.array(o.f.line_difference[0]).std()}")
    plot_all(o)
    plot_distance(o)
    # plot_sigmas(o)
