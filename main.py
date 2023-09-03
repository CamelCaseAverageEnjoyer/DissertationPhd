from main_objects import *
from plot_func import *


if __name__ == '__main__':
    # Инициализация
    start_navigation_tolerance = 0.9
    t_integrate = 1e4
    dt = 10.

    start_navigation = ['perfect', 'near', 'random'][0]
    model_c = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
    method_navigation = ['kalman_filter rv', 'kalman_filter rvq', 'kalman_filter rv all'][0]
    kalman_coef = [[0.0001, 0.01, 1e-09, 1e-12, 1e-09, 1e-08], ['d', 'q', 'p1', 'p2', 'p3', 'r']]  # посчитанные
    # kalman_coef = [[1e-4, 3e-7, 1e-15, 1e-15, 1e-15, 1e-2], ['d', 'q', 'p1', 'p2', 'p3', 'r']]  # 1
    # kalman_coef = [[1e-3, 1e-15, 1e-6, 1e-10, 1e-10, 1e-2], ['d', 'q', 'p1', 'p2', 'p3', 'r']]  # N

    o = Objects(n_c=1, n_f=1, kalman_coef=kalman_coef,
                model_c=model_c, dt=dt, start_navigation=start_navigation,
                start_navigation_tolerance=start_navigation_tolerance, method_navigation=method_navigation)
    o.f.gain_mode = ['isotropic', 'ellipsoid'][1]
    o.f.ellipsoidal_signal = 0.99
    o.p.show_rate = 1  # int(round(1 / dt))
    o.p.kalman_single_search = True
    o.p.is_aero = False

    # Интегрирование
    o.p.integrate(t=t_integrate)

    # Результаты
    plot_all(o)
    plot_distance(o)
    # plot_signals(o)
