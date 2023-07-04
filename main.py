from main_objects import *
from plot_func import *


if __name__ == '__main__':
    # Инициализация
    o = Objects(n_c=1, n_f=5, model_c='1U', dt=10.)
    o.p.show_rate = 1
    o.f.ellipsoidal_signal = 0.9

    # Интегрирование + результаты
    o.p.integrate(t=1e4)
    plot_distance(o)
    # plot_signals(o)
    # plot_all(o)
