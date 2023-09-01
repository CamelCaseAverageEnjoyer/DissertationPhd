from main_objects import *
from plot_func import *


if __name__ == '__main__':
    # Инициализация
    start_navigation_tolerance = 0.9
    t_integrate = 5e3
    dt = 1.

    start_navigation = ['perfect', 'near', 'random'][1]
    model_c = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
    method_navigation = ['kalman_filter rv', 'kalman_filter rvq', 'kalman_filter rv all'][2]
    kalman_coef = [[1e-3, 1e-12, 1e-6, 1e-10, 1e-10, 1e-2],
                   ['d', 'q', 'p1', 'p2', 'p3', 'r']]

    o = Objects(n_c=1, n_f=2, kalman_coef=kalman_coef,
                model_c=model_c, dt=dt, start_navigation=start_navigation,
                start_navigation_tolerance=start_navigation_tolerance, method_navigation=method_navigation)
    o.f.gain_mode = ['isotropic', 'ellipsoid'][1]
    o.f.ellipsoidal_signal = 0.99
    o.p.show_rate = 10  # int(round(1 / dt))
    o.p.kalman_single_search = True
    o.p.is_aero = False

    # Интегрирование
    o.p.integrate(t=t_integrate)

    # Результаты
    plot_all(o)
    # plot_distance(o)
    # plot_signals(o)
