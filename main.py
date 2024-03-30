from main_objects import *


if __name__ == '__main__':
    # Инициализация
    start_navigation_tolerance = 0.1
    t_integrate = 1e5
    dt = 10.

    start_navigation = ['perfect', 'near', 'random'][1]
    model_c = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
    method_navigation = ['kalman_filter rv', 'kalman_filter rvq', 'kalman_filter rv all', 'kalman_filter rvq all'][2]
    kalman_coef = {'q': [1e-10, 1e-8], 'p': [1e-8] * 4, 'r': 1e1}

    o = Objects(n_c=1, n_f=10, kalman_coef=kalman_coef,
                model_c=model_c, dt=dt, start_navigation=start_navigation,
                start_navigation_tolerance=start_navigation_tolerance, method_navigation=method_navigation)
    o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas'][0]
    o.c.gain_mode = o.f.gain_mode
    o.p.k.gain_mode = o.f.gain_mode
    o.p.is_aero = False
    o.p.show_rate = 1

    # Интегрирование
    o.p.integrate(t=t_integrate, animate=False)

    # Результаты
    tmp = np.array([np.linalg.norm(o.f.line_difference[0][i]) for i in range(len(o.f.line_difference[0]))])
    print(f"Математическое ожидание: {tmp.mean()}, Среднее отклонение: {tmp.std()}")
    # plot_all(o)
    plot_distance(o)
    # plot_sigmas(o)
