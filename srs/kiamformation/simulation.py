"""Здесь много всего так-то"""
import matplotlib.pyplot as plt
import time

from config import *

def timer(func):
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f'Выполнено "{func.__name__}" за {run_time:.4f} секунд')
        return value
    return wrapper_timer

def save_simulation_trajectories(o: Objects, text: str):
    o.p.record.to_csv(f'{text}.csv', index=False, sep=";")
    my_print(f"В файл {text} записаны траектории", color='y')

def load_simulation_trajectories(o: Objects, text: str):
    from pandas import read_csv
    o.p.record = read_csv(text, sep=";")
    my_print(f"Из файла {text} прочитаны траектории", color='y')

def find_mirror_solutions(o: Objects):
    """
    Пока что сделано только для изотропных антенн, без учёта ориентации
    :param o:
    :return:
    """
    from dynamics import get_rand_c
    # c_0 = o.f.c_hkw
    c_r = get_rand_c(v=o.v)
    pass  # o.v.MEASURES_VECTOR

def solve_minimization_new(o: Objects, config_choose_n: int):
    from scipy.optimize import minimize
    global o_rand, measure_lst
    t = o.p.record['t'].to_numpy()
    n_measure = int(o.p.record.loc[1, 'MEASURES_VECTOR N'])
    measures_real = [o.p.record[f'MEASURES_VECTOR {i}'].to_numpy() for i in range(n_measure)]

    # Инициализация
    o.v.IF_ANY_PRINT = o.v.IF_TEST_PRINT = False
    o_rand = Objects(v=o.v)
    o_rand.c = o.c
    # Подгон
    for i in range(o_rand.f.n):
        o_rand.f.r_orf[i] = o.f.apriori_params['r orf'][i]  # o.f.r_orf[i].copy()
        o_rand.f.v_orf[i] = o.f.apriori_params['v orf'][i]  # o.f.v_orf[i].copy()
    o_rand.f.update_c(v=o.v)
    # Запись приближения
    if np.any(list(o.v.DYNAMIC_MODEL.values())):  # Если J2 или aero drag
        x0 = np.append(o_rand.f.r_orf[0].flatten(), o_rand.f.v_orf[0].flatten())
    else:
        x0 = o_rand.f.c_hkw[0]

    def get_error(x):
        """Функция для вызова scipy.optimize.minimize"""
        global o_rand, measure_lst
        o_rand.reset(config_choose_n)
        o_rand.c = o.c
        o.v.IF_NAVIGATION = False

        if np.any(list(o.v.DYNAMIC_MODEL.values())):  # Если J2 или aero drag
            o_rand.f.r_orf[0] = x[:len(x)//2].copy()
            o_rand.f.v_orf[0] = x[len(x)//2:].copy()
            o_rand.f.param_fitting(v=o_rand.v)
            o_rand.f.update_c(v=o_rand.v)
        else:
            o_rand.f.c_hkw[0] = x.copy()
        o_rand.integrate(t=o.p.t)
        # print(f"o2: r={o_rand.f.r_orf[0]}, v={o_rand.f.v_orf[0]}")

        measures = [o_rand.p.record[f'MEASURES_VECTOR {i}'].to_numpy() for i in range(n_measure)]
        for i in range(n_measure):
            measure_lst[i].append(measures[i])

        return sum([np.abs(measures[i] - measures_real[i]).sum() for i in range(n_measure)])

    measure_lst = [[] for _ in range(len(measures_real))]
    # res = minimize(get_error, x0, method='trust-constr', options={'verbose': 3})
    # res = minimize(get_error, x0, method='COBYLA', options={'tol': 1e-15, 'disp': True})
    res = minimize(get_error, x0, method='TNC', options={'ftol': 1e-10, 'disp': True})
    print(f"get:  {res.x}")
    print(f"need: {o.f.c_hkw[0]}")
    print(f"----------------\nВсе результаты: {res}")

    # Отображение - позже занести measures_est, r_est в цикл, считать для каждого чипсата
    fig, axs = plt.subplots(1, 2, figsize=(20, 8))

    for i in range(n_measure):
        for im, m in enumerate(measure_lst):
            for ik, yk in enumerate(m):
                axs[0].plot(t, yk, ":", label='scipy history' if ik + im == 0 else None, c='palegreen')
        y1 = measures_real[i]
        y2 = o_rand.p.record[f'MEASURES_VECTOR {i}'].to_numpy()
        axs[0].plot(t, y1, "-", label='real', c=o.v.MY_COLORS[i])
        axs[0].plot(t, y2, "-", label='estimate', c=o.v.MY_COLORS[i+1])
        axs[0].plot(t, np.abs(y1 - y2), "-", label='error', c='k')


    for i, c in enumerate('xyz'):
        for k in range(o.f.n):
            y1 = o.p.record[f'{o.f.name} r {c} orf {k}'].to_numpy()
            y2 = o_rand.p.record[f'{o.f.name} r {c} orf {k}'].to_numpy()
            axs[1].plot(t, y1, "-", label=f'(id={k}) {c} real', c=o.v.MY_COLORS[i])
            axs[1].plot(t, y2, ":", label=f'(id={k}) {c} estimate', c=o.v.MY_COLORS[i])

    for j in range(2):
        axs[j].set_title(["Измерение", "Состояние"][j])
        axs[j].grid(True)
        axs[j].legend()
    plt.show()
