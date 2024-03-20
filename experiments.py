"""Здесь функции расчёта красивых графиков для красивых презентаций"""
from typing import List, Any

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
from main_objects import *
from tiny_functions import *

def do_tolerance_table(t_integrate: float = 1e4, dt: float = 10., aero: bool = False, repeat: int = 3) -> None:
    """Богом гонимая таблца точностей при разном количестве аппаратов, при разных точностях
    :argument:
    - t_integrate -- время интегрирования уравнений СИ (default 1e4);               \n
    - dt -- шаг по времени СИ (default 1e1);                                        \n
    - aero -- флаг учёта аэродинамического сопротивления (default False);           \n
    - repeat -- количество генераций движения на одну ячейку таблицы (default 3);   \n
    :return: None
    """
    kalman_coef = {'q': 1e-15, 'p': [1e-8] * 3, 'r': 1e0}
    tolerance_list = [0.8, 0.9, 1.]
    n_f_list = [10, 20, 30]
    counter = 0
    time_begin = datetime.now()
    n_total = len(tolerance_list) * len(tolerance_list) * repeat
    res_mean = [[0. for _ in range(len(n_f_list))] for _ in range(len(tolerance_list))]
    res_std = [[0. for _ in range(len(n_f_list))] for _ in range(len(tolerance_list))]

    for i_t in range(len(tolerance_list)):
        for i_n in range(len(n_f_list)):
            for _ in range(repeat):
                # Отображение процесса
                counter += 1
                per = int(math.ceil(30 * counter / n_total))
                real_per = int(math.ceil(100 * counter / n_total))
                print(f"{real_per}% [{'#' * per + ' ' * (30 - per)}], время: {datetime.now() - time_begin}")

                # Сам процесс
                o = Objects(n_c=1, n_f=n_f_list[i_n], kalman_coef=kalman_coef,
                            model_c='1U', dt=dt, start_navigation='near', if_any_print=False,
                            start_navigation_tolerance=tolerance_list[i_t], method_navigation='kalman_filter rvq all')
                o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas'][3]
                o.c.gain_mode = o.f.gain_mode
                o.p.k.gain_mode = o.f.gain_mode
                o.p.is_aero = aero
                o.p.integrate(t=t_integrate)

                tmp_mean = 0.
                tmp_std = 0.
                for i in range(n_f_list[i_n]):
                    tmp_mean += np.array(o.f.line_difference[i]).mean()
                    tmp_std += np.array(o.f.line_difference[i]).std()
                res_mean[i_t][i_n] += tmp_mean / repeat / n_f_list[i_n]
                res_std[i_t][i_n] += tmp_std / repeat / n_f_list[i_n]

    tmp_col = {i: f"{int(tolerance_list[i] * 100)}%" for i in range(len(tolerance_list))}
    tmp_row = {i: n_f_list[i] for i in range(len(n_f_list))}
    res_mean_df = pd.DataFrame(res_mean).rename(columns=tmp_col, index=tmp_row)
    res_std_df = pd.DataFrame(res_std).rename(columns=tmp_col, index=tmp_row)
    print(f"Средние средние: \n{res_mean_df}\nСредние отклонения: \n{res_std_df}")

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
    sns.heatmap(res_mean_df, annot=True, fmt=".1f", ax=axs[0], cmap=cmap)
    sns.heatmap(res_std_df, annot=True, fmt=".2f", ax=axs[1], cmap=cmap)
    axs[0].set_title("Погрешность определения положения")
    axs[1].set_title("Отклонение погрешности навигации")
    for i in range(2):
        axs[i].set_xlabel("Точность приорной информации")
        axs[i].set_ylabel("Количество чипсатов")
    plt.show()

def params_search(method_navigation: str = 'kalman_filter rvq', t_integrate: float = 1e4, dt: float = 10.,
                  aero: bool = False) -> list:
    """Перебирает параметры фильтра Калмана в святой надежде на светлое будущее
    :argument:
    - method_navigation -- метод наблюдения (default 'kalman_filter rvq');  \n
    - t_integrate -- время интегрирования уравнений СИ (default 1e4);       \n
    - dt -- шаг по времени СИ (default 1e1);                                \n
    - aero -- флаг учёта аэродинамического сопротивления (default False);   \n
    :return: best_params -- [q, p1, p2, p2, r]"""
    min_discrepancy = 1e10
    count = 0
    best_params = []
    q_list = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    r_list = [1e-4, 1e-2]
    p1_list = [1e-12, 1e-10, 1e-8]
    p2_list = [1e-12, 1e-10, 1e-8]
    N = len(q_list) * len(r_list) * len(p1_list) * len(p2_list)
    time_begin = datetime.now()
    for q in q_list:
        for r in r_list:
            for p1 in p1_list:
                for p2 in p2_list:
                    tmp = 0
                    count += 1
                    if count % 10 == 0:
                        print(f"Поиск коэффициентов фильтра Калмана {count}/{N}={int(100 * count / N)}%, " +
                              real_workload_time(n=count, n_total=N,
                                                 time_begin=time_begin, time_now=datetime.now()))
                    for _ in range(1):
                        o = Objects(n_c=1, n_f=1,
                                    kalman_coef={'q': q, 'p': [p1, p2, p2], 'r': r},
                                    model_c='1U', dt=dt, start_navigation='near',  # random perfect
                                    start_navigation_tolerance=0.98,
                                    method_navigation=method_navigation, if_any_print=False)
                        o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas'][3]
                        o.f.ellipsoidal_signal = 0.99
                        o.p.is_aero = aero
                        o.p.integrate(t=t_integrate)

                        tmp += np.sum(np.abs(np.array(o.c.real_dist[0][0]) - np.array(o.c.calc_dist[0][0])))
                    if tmp < min_discrepancy:
                        min_discrepancy = tmp
                        best_params = [q, p1, p2, p2, r]
    print(f"Подходящие параметры: 'q': {best_params[0]}, 'p': {best_params[1:4]}, 'r': {best_params[4]}")
    print(f"best_params={best_params}")
    return best_params

def solve_minimization(aero: bool = False, quaternion: bool = False, quaternion_but_i_dont_give_a_fuck: bool = False,
                       random_rate: float = 0.1, post_factor: int = 1, n: int = 1, r_noise: float = 0,
                       dt: float = 1, t_between: float = 20, n_points: int = 5):
    """Функция берёт и суёт вашу обратную задачу в Scipy"""
    import scipy

    Radius_orbit = 6800e3
    mu = 5.972e24 * 6.67408e-11
    w = np.sqrt(mu / Radius_orbit ** 3)
    v_orb = np.sqrt(mu / Radius_orbit)
    h_orb = Radius_orbit - 6371e3

    # Задание параметров
    show_factor = 1000
    t = [60 * t_between * i for i in range(n_points)]
    colors = ["violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
              "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
              "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
              "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
              "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
              "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
              "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
              "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
              "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
              "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine"]

    # Генерация движения
    def get_rand_rvs(r_spread: float = 100, v_spread: float = 0.01):
        r_lcl = np.random.uniform(-r_spread, r_spread, 3)
        v_lcl = np.random.uniform(-v_spread, v_spread, 3)
        c_lcl = np.random.uniform(0, 0.05 ** 2, 1)
        return np.append(np.append(r_lcl, v_lcl), c_lcl).tolist()

    rvs = []
    rvs_rand = []
    C = []
    C_real = []
    C_rand = []
    ev = np.array([1., 0., 0.])
    len_x = 9 if quaternion or quaternion_but_i_dont_give_a_fuck else 6
    len_x_res = 9 if quaternion else 6
    if aero:
        for i in range(n):
            rvs += get_rand_rvs()
            rvs_rand += get_rand_rvs()
    else:
        for i in range(n):
            tmp = get_rand_c(w=w, if_quaternion=quaternion or quaternion_but_i_dont_give_a_fuck)
            C += tmp
            C_real += tmp[0:len_x_res]
            C_rand += get_rand_c(if_no=quaternion_but_i_dont_give_a_fuck and not quaternion, w=w,
                                 if_quaternion=quaternion or quaternion_but_i_dont_give_a_fuck)

    # Расчёт измерений
    R = []
    if aero:
        rvs_copy = rvs.copy()
        for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
            r = []
            time = j * dt
            for i in range(n):
                rvs_copy[i * 7:(i + 1) * 7] = get_integrate(rvs=rvs_copy[i * 7:(i + 1) * 7], dt=dt, h_orb=h_orb,
                                                            v_orb=v_orb, w=w)
            for t1 in t:
                if time == t1:
                    for i in range(n):
                        r += rvs_copy[i * 7:i * 7 + 3]
                    for i1 in range(n):
                        for i2 in range(i1):
                            dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                            R += [dist_real * (1 + np.random.uniform(-r_noise, r_noise))]
    else:
        for j in t:
            r = []
            for i in range(n):
                r += [np.array(r_hkw(C[i * len_x_res:i * len_x_res + 6], w, j))]
            for i1 in range(n):
                for i2 in range(i1 + 1):
                    dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                    if quaternion or quaternion_but_i_dont_give_a_fuck:
                        g1 = get_gain(euler2rot_matrix(C[i1 * 9 + 6] / 30, C[i1 * 9 + 7] / 30, C[i1 * 9 + 8] / 30) @ ev)
                        g2 = get_gain(euler2rot_matrix(C[i2 * 9 + 6] / 30, C[i2 * 9 + 7] / 30, C[i2 * 9 + 8] / 30) @ ev)
                        signal_rate = g1 if i1 == i2 else g1 * g2

                        dist_calc = dist_real * np.sqrt(signal_rate)
                    else:
                        dist_calc = dist_real
                    R += [dist_calc * (1 + np.random.uniform(-r_noise, r_noise))]

    # Отображение
    x, y, z = ([], [], [])
    rvs_copy = rvs.copy()
    ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')

    if aero:
        for i in range(n):
            x_, y_, z_ = ([], [], [])
            for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
                time = j * dt
                rvs_copy[i * 7:(i + 1) * 7] = get_integrate(rvs_copy[i * 7:(i + 1) * 7])
                x_ += [rvs_copy[i * 7 + 0]]
                y_ += [rvs_copy[i * 7 + 1]]
                z_ += [rvs_copy[i * 7 + 2]]
                for t1 in t:
                    if time == t1:
                        x += [rvs_copy[i * 7 + 0]]
                        y += [rvs_copy[i * 7 + 1]]
                        z += [rvs_copy[i * 7 + 2]]
            ax.plot(x_, y_, z_, c=colors[i])
    else:
        for i in range(n):
            for j in t:
                tmp = r_hkw(C[i * len_x:i * len_x + 6], w, j)
                x += [tmp[0]]
                y += [tmp[1]]
                z += [tmp[2]]
            x_, y_, z_ = ([], [], [])
            for j in np.linspace(t[0], t[-1], 100):
                tmp = r_hkw(C[i * len_x:i * len_x + 6], w, j)
                x_ += [tmp[0]]
                y_ += [tmp[1]]
                z_ += [tmp[2]]
            ax.plot(x_, y_, z_, c=colors[i])
    ax.scatter(x, y, z)
    plt.show()

    def local_func(arg: any):
        loss = 0
        counter = 0
        if aero:
            rvs_copy = arg.copy().tolist()
            for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
                r = []
                time = j * dt
                for i in range(n):
                    rvs_copy[i * 7:(i + 1) * 7] = get_integrate(rvs_copy[i * 7:(i + 1) * 7])
                for t1 in t:
                    if time == t1:
                        for i in range(n):
                            r += rvs_copy[i * 7:i * 7 + 3]
                        for i1 in range(n):
                            for i2 in range(i1):
                                R_ = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                                loss += (R_ - R[counter]) ** 2
                                counter += 1
        else:
            C_ = arg
            for j in t:
                r = []
                for i in range(n):
                    r += [np.array(r_hkw(C_[i * len_x_res:i * len_x_res + 6], w, j))]
                for i1 in range(n):
                    for i2 in range(i1 + 1):
                        dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                        if quaternion and not quaternion_but_i_dont_give_a_fuck:
                            g1 = get_gain(
                                euler2rot_matrix(C_[i1 * 9 + 6] / 30, C_[i1 * 9 + 7] / 30, C_[i1 * 9 + 8] / 30) @ ev)
                            g2 = get_gain(
                                euler2rot_matrix(C_[i2 * 9 + 6] / 30, C_[i2 * 9 + 7] / 30, C_[i2 * 9 + 8] / 30) @ ev)
                            signal_rate = g1 if i1 == i2 else g1 * g2
                            dist_calc = dist_real * np.sqrt(signal_rate)
                        else:
                            dist_calc = dist_real
                        R_ = dist_calc * (1 + np.random.uniform(-r_noise, r_noise))

                        loss += (R_ - R[counter]) ** 2
                        counter += 1
        return loss / counter

    # Расчёт
    if aero:
        x0 = np.array(rvs_rand) * random_rate + np.array(rvs) * (1 - random_rate)
    else:
        x0 = np.array(C_rand) * random_rate + np.array(C_real) * (1 - random_rate)
    res = scipy.optimize.minimize(local_func, x0, tol=1e-4, method='trust-constr', options={'verbose': 3})
    print(f"Результатики: {res.x}\nА надо бы {C}\n")

    # Повторное отображение
    if aero:
        ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
        x, y, z = ([], [], [])
        rvs_copy = rvs.copy()
        rvs_res = res.x.copy()
        a = [[] for i in range(n)]
        for i in range(n):
            x_, y_, z_ = ([], [], [])
            x_1, y_1, z_1 = ([], [], [])
            for j in np.linspace(0, t[-1] * post_factor, show_factor * int(t[-1] * post_factor // dt), endpoint=False):
                time = j * dt
                rvs_res[i * 7:(i + 1) * 7] = get_integrate(rvs_res[i * 7:(i + 1) * 7])
                x_ += [rvs_res[i * 7 + 0]]
                y_ += [rvs_res[i * 7 + 1]]
                z_ += [rvs_res[i * 7 + 2]]
                rvs_copy[i * 7:(i + 1) * 7] = get_integrate(rvs_copy[i * 7:(i + 1) * 7])
                x_1 += [rvs_copy[i * 7 + 0]]
                y_1 += [rvs_copy[i * 7 + 1]]
                z_1 += [rvs_copy[i * 7 + 2]]
                a[i] += [np.linalg.norm(np.array(rvs_res[i * 7:i * 7 + 3]) - np.array(rvs_copy[i * 7:i * 7 + 3]))]
                for t1 in t:
                    if time == t1:
                        x += [rvs_copy[i * 7 + 0]]
                        y += [rvs_copy[i * 7 + 1]]
                        z += [rvs_copy[i * 7 + 2]]
            ax.plot(x_, y_, z_, c=colors[i])
            ax.plot(x_1, y_1, z_1, c=colors[i], ls=":")
        ax.scatter(x, y, z)
        ax.set_title("Траектории дейтвительные и расчитанные")
        ax.set_xlabel("x, м")
        ax.set_ylabel("y, м")
        ax.set_zlabel("z, м")
        plt.show()

        for i in range(n):
            plt.plot(np.linspace(0, t[-1] * post_factor, show_factor * int(t[-1] * post_factor // dt), endpoint=False),
                     a[i], c=colors[i])
        plt.title("Реальная ошибка")
        plt.xlabel("Время, с")
        plt.ylabel("Ошибка, м")
        plt.show()

        # Расчёт и отображение измерений
        loss = []
        counter = 0
        rvs_res = res.copy()
        for j in np.linspace(0, t[-1] * post_factor, int(t[-1] * post_factor // dt), endpoint=False):
            tmp = 0.
            r = []
            time = j * dt
            for i in range(n):
                rvs_res[i * 7:(i + 1) * 7] = get_integrate(rvs_res[i * 7:(i + 1) * 7])
            for t1 in t[0:-1]:
                if time == t1:
                    for i in range(n):
                        r += rvs_res[i * 7:i * 7 + 3]
                    for i1 in range(n):
                        for i2 in range(i1):
                            # print(f"{counter+1} / {len(R)}, t = {j}")
                            R_ = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                            tmp += (R_ - R[counter]) ** 2
                            counter += 1
            loss += [np.sqrt(tmp) / counter]
        for j in t:
            plt.plot([j, j], [np.min(loss), np.max(loss)], c='gray')
        plt.plot([t[0], t[-1] * post_factor], [0, 0], c='gray')
        plt.plot(np.linspace(0, t[-1] * post_factor, int(t[-1] * post_factor // dt), endpoint=False), loss)
        plt.title("Расчётная ошибка")
        plt.xlabel("Время, с")
        plt.ylabel("Ошибка, м")
        plt.show()
    else:
        ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
        # x, y, z = ([], [], [])
        for i in range(n):
            '''for j in t:
                tmp = r_hkw(C[i*len_x:i*len_x+6], w, j)
                x += [tmp[0]]
                y += [tmp[1]]
                z += [tmp[2]]'''
            x_, y_, z_ = ([], [], [])
            for j in np.linspace(t[0], t[-1] * post_factor, show_factor):
                tmp = r_hkw(C[i * len_x:i * len_x + 6], w, j)
                x_ += [tmp[0]]
                y_ += [tmp[1]]
                z_ += [tmp[2]]
            ax.plot(x_, y_, z_, c=colors[i])
            x_, y_, z_ = ([], [], [])
            for j in np.linspace(t[0], t[-1] * post_factor, show_factor):
                tmp = r_hkw(res.x[i * len_x_res:i * len_x_res + 6], w, j)
                x_ += [tmp[0]]
                y_ += [tmp[1]]
                z_ += [tmp[2]]
            ax.plot(x_, y_, z_, c=colors[i], ls=":")
        # ax.scatter(x, y, z)
        ax.set_title("Траектории дейтвительные и расчитанные")
        ax.set_xlabel("x, м")
        ax.set_ylabel("y, м")
        ax.set_zlabel("z, м")
        plt.show()

        max_value = 0.
        for i in range(n):
            a = []
            for j in np.linspace(t[0], t[-1] * post_factor, 100):
                tmp = np.array(r_hkw(C[i * len_x:i * len_x + 6], w, j)) - np.array(
                    r_hkw(res.x[i * len_x_res:i * len_x_res + 6], w, j))
                a += [np.linalg.norm(tmp)]
            plt.plot(np.linspace(t[0], t[-1] * post_factor, 100), a, c=colors[i])
            max_value = max(max_value, np.max(a))
        plt.title("Реальная ошибка")
        plt.xlabel("Время, с")
        plt.ylabel("Ошибка, м")
        plt.show()

        # Расчёт и отображение измерений ТУТ НЕТ УЧЁТА ПОВОРОТА и ДИАГРАММ НАПРАВЛЕННОСТЕЙ
        loss = []
        C_ = res.x.tolist()
        for j in np.linspace(t[0], t[-1] * post_factor, 200):
            counter = 0
            tmp = 0
            r = []
            r_ = []
            for i in range(n):
                r += [np.array(r_hkw(C[i * len_x:i * len_x + 6], w, j))]
                r_ += [np.array(r_hkw(res.x[i * len_x_res:i * len_x_res + 6], w, j))]
            for i1 in range(n):
                for i2 in range(i1):
                    dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                    if quaternion or quaternion_but_i_dont_give_a_fuck:
                        g1 = get_gain(euler2rot_matrix(C[i1 * 9 + 6] / 30, C[i1 * 9 + 7] / 30, C[i1 * 9 + 8] / 30) @ ev)
                        g2 = get_gain(euler2rot_matrix(C[i2 * 9 + 6] / 30, C[i2 * 9 + 7] / 30, C[i2 * 9 + 8] / 30) @ ev)
                        signal_rate = g1 if i1 == i2 else g1 * g2

                        dist_calc = dist_real * np.sqrt(signal_rate)
                    else:
                        dist_calc = dist_real
                    R = dist_calc

                    dist_real = np.linalg.norm(r_[i1]) if i1 == i2 else np.linalg.norm(r_[i1] - r_[i2])
                    if quaternion and not quaternion_but_i_dont_give_a_fuck:
                        g1 = get_gain(euler2rot_matrix(C_[i1 * 9 + 6] / 30, C_[i1 * 9 + 7] / 30, C_[i1 * 9 + 8] / 30) @ ev)
                        g2 = get_gain(euler2rot_matrix(C_[i2 * 9 + 6] / 30, C_[i2 * 9 + 7] / 30, C_[i2 * 9 + 8] / 30) @ ev)
                        signal_rate = g1 if i1 == i2 else g1 * g2
                        dist_calc = dist_real * np.sqrt(signal_rate)
                    else:
                        dist_calc = dist_real
                    R_ = dist_calc

                    tmp += (R_ - R) ** 2
                    counter += 1
            loss += [np.sqrt(tmp / counter)]
        for j in t:
            plt.plot([j, j], [0, max_value], c='gray')
        plt.plot([t[0], t[-1] * post_factor], [0, 0], c='gray')
        plt.plot(np.linspace(t[0], t[-1] * post_factor, 200), loss)
        plt.title("Расчётная ошибка")
        plt.xlabel("Время, с")
        plt.ylabel("Ошибка, м")
        plt.show()

        for j in t:
            plt.plot([j, j], [0, np.max(loss)], c='gray')
        plt.plot([t[0], t[-1] * post_factor], [0, 0], c='gray')
        plt.plot(np.linspace(t[0], t[-1] * post_factor, 200), loss)
        plt.title("Расчётная ошибка")
        plt.xlabel("Время, с")
        plt.ylabel("Ошибка, м")
        plt.show()


if __name__ == '__main__':
    pass
    # params_search()
    # do_tolerance_table()
