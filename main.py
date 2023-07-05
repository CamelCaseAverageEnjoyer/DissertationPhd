from main_objects import *
from plot_func import *


if __name__ == '__main__':
    # Инициализация
    t_integrate = 2e5  # e5
    dt = 1.

    o = Objects(n_c=1, n_f=1, model_c='1U', dt=dt)
    o.f.ellipsoidal_signal = 0.99
    o.p.show_rate = 10  # int(round(1 / dt))
    o.p.kalman_single_search = True
    o.p.is_aero = True

    # Интегрирование + результаты
    o.p.integrate(t=t_integrate)
    plot_distance(o)
    plot_all(o)
    # plot_signals(o)
