from main_objects import *
from plot_func import *


if __name__ == '__main__':
    o = Objects(n_c=1, n_f=5, model_c='1U', dt=10., show_rate=1)
    o.p.integrate(1e5)
    plot_signals(o)
    # plot_all(o)
