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
