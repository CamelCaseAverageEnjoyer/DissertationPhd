from colorama import Fore, Style
from typing import Union
from numpy.random import uniform
from datetime import datetime
# import scipy
# import numpy as np
# import kiam_astro

from tiny_functions import *


class Cosmetic:
    def __init__(self, if_any_print: bool = False):
        self.if_any_print = if_any_print

    def my_print(self, txt, color=None) -> None:
        if self.if_any_print:
            if color is None:
                print(Style.RESET_ALL + txt)
            if color == "b":
                print(Fore.BLUE + txt + Style.RESET_ALL)
            if color == "g":
                print(Fore.GREEN + txt + Style.RESET_ALL)
            if color == "y":
                print(Fore.YELLOW + txt + Style.RESET_ALL)
            if color == "r":
                print(Fore.RED + txt + Style.RESET_ALL)
            if color == "c":
                print(Fore.CYAN + txt + Style.RESET_ALL)
            if color == "m":
                print(Fore.MAGENTA + txt + Style.RESET_ALL)

class FemtoSat:
    def __init__(self, w_orb: float, n: int = 10, r_spread: float = 1e2, v_spread: float = 1e-1,
                 start_navigation: str = 'near', start_navigation_tolerance: float = 0.9, w_spread: float = 1e-3):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        # Предопределённые параметры
        self.gain_modes = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
        self.navigations = ['perfect', 'near', 'random']

        # Общие параметры
        self.n = n
        self.size = [0.03, 0.03]
        self.mass = 0.03
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001

        self.ellipsoidal_signal = 0.9
        self.gain_mode = self.gain_modes[1]

        # Индивидуальные параметры
        self.r_orf = [uniform(-r_spread, r_spread, 3) for _ in range(self.n)]
        self.v_orf = [uniform(-v_spread, v_spread, 3) for _ in range(self.n)]
        self.w_orf = [uniform(-1., 1., 3) for _ in range(self.n)]
        self.q = [uniform(-1, 1, 4) for _ in range(self.n)]
        self.q_ = [uniform(-1, 1, 4) for _ in range(self.n)]
        for i in range(self.n):
            self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
            self.w_orf[i] *= uniform(-w_spread, w_spread) / np.linalg.norm(self.w_orf[i])
        self.c_hkw = [[] for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
            self.q_[i] /= np.linalg.norm(self.q_[i])
        self.line = [[] for _ in range(self.n)]
        self.line_kalman = [[] for _ in range(self.n)]
        self.line_difference = [[] for _ in range(self.n)]
        self.z_difference = [[] for _ in range(self.n)]
        self.signal_power = [[[] for _ in range(self.n)] for _ in range(self.n)]
        self.real_dist = [[[] for _ in range(self.n)] for _ in range(self.n)]
        self.calc_dist = [[[] for _ in range(self.n)] for _ in range(self.n)]
        self.calc_dist_ = [[[] for _ in range(self.n)] for _ in range(self.n)]
        tmp_poor = [np.array([uniform(-r_spread, r_spread) for _ in range(3)] +
                             [uniform(-v_spread, v_spread) for _ in range(3)] +
                             list(self.q_[i][1:4])) for i in range(self.n)]
        tmp_good = [np.append(np.append(self.r_orf[i], self.v_orf[i]), list(self.q[i][1:4])) for i in range(self.n)]

        # Параметры начального приближения
        if start_navigation == self.navigations[0]:
            self.rv_orf_calc = tmp_good
        elif start_navigation == self.navigations[1]:
            self.rv_orf_calc = [tmp_good[i] * start_navigation_tolerance +
                                tmp_poor[i] * (1 - start_navigation_tolerance) for i in range(self.n)]
        else:
            self.rv_orf_calc = tmp_poor

    def get_gain(self, r: Union[float, np.ndarray], mode3: bool = False):
        r1 = r / np.linalg.norm(r)
        if self.gain_mode == self.gain_modes[1]:
            return np.linalg.norm([r1[0] * 1, r1[1] * self.ellipsoidal_signal, r1[2] * self.ellipsoidal_signal])
        if self.gain_mode == self.gain_modes[2]:
            return 1 - np.dot(r1, np.array([0, 0, 1]))**2
        if self.gain_mode == self.gain_modes[3] or (mode3 and self.gain_mode == self.gain_modes[4]):
            return 2 - np.dot(r1, np.array([1, 0, 0]))**2 - np.dot(r1, np.array([0, 1, 0]))**2
        if self.gain_mode == self.gain_modes[4] and not mode3:
            return [1 - np.dot(r1, np.array([1, 0, 0]))**2, 1 - np.dot(r1, np.array([0, 1, 0]))**2]
        return 1

class CubeSat:
    def __init__(self, w_orb: float, n_f: int, n: int = 1, model: str = '1U', r_spread: float = 1e-1,
                 v_spread: float = 1e-4):
        """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.\n
        Все величны представлены в СИ."""
        # Предопределённые параметры
        self.gain_modes = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
        models = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        masses = [2., 3., 4., 6., 12., 24.]
        mass_center_errors = [[0.02, 0.02, 0.02], [0.02, 0.02, 0.03], [0.02, 0.02, 0.045],
                              [0.02, 0.02, 0.07], [4.5, 2., 7.], [4.5, 4.5, 7.]]
        sizes = [[0.1, 0.1, 0.1135], [0.1, 0.1, 0.1702], [0.1, 0.1, 0.227],
                 [0.1, 0.1, 0.3405], [0.2263, 0.1, 0.366], [0.2263, 0.2263, 0.366]]

        # Общие параметры
        self.n = n
        self.model = model
        self.model_number = models.index(model)
        self.mass = masses[self.model_number]
        self.mass_center_error = mass_center_errors[self.model_number]
        self.r_mass_center = np.array([uniform(-i, i) for i in self.mass_center_error])
        self.size = sizes[self.model_number]

        self.ellipsoidal_signal = 0.9
        self.gain_mode = self.gain_modes[0]

        # Индивидуальные параметры
        self.r_orf = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.v_orf = [np.array([uniform(-v_spread, v_spread) for _ in range(3)]) for _ in range(self.n)]
        for i in range(self.n):
            self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
        self.c_hkw = [[] for _ in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]
        self.signal_power = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.real_dist = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.calc_dist = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.calc_dist_ = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.kalm_dist = [[[] for _ in range(n_f)] for _ in range(self.n)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

    def get_gain(self, r: Union[float, np.ndarray], mode3: bool = False):
        r1 = r / np.linalg.norm(r)
        if self.gain_mode == self.gain_modes[1]:
            return np.linalg.norm([r1[0] * 1, r1[1] * self.ellipsoidal_signal, r1[2] * self.ellipsoidal_signal])
        if self.gain_mode == self.gain_modes[2]:
            return 1 - np.dot(r1, np.array([0, 0, 1]))**2
        if self.gain_mode == self.gain_modes[3] or (mode3 and self.gain_mode == self.gain_modes[4]):
            return 2 - np.dot(r1, np.array([1, 0, 0]))**2 - np.dot(r1, np.array([0, 1, 0]))**2
        if self.gain_mode == self.gain_modes[4] and not mode3:
            return [1 - np.dot(r1, np.array([1, 0, 0]))**2, 1 - np.dot(r1, np.array([0, 1, 0]))**2]
        return 1

class KalmanFilter:
    def __init__(self, f: FemtoSat, c: CubeSat, s: Cosmetic, p, dt: float, w_orb: float, kalman_coef: dict,
                 orientation: bool = False, single_femto_filter: bool = True):
        self.gain_modes = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
        self.gain_mode = self.gain_modes[0]
        self.f = f
        self.c = c
        self.s = s
        self.p = p
        self.dt = dt
        self.q_coef = kalman_coef['q']
        self.p_coef = kalman_coef['p']
        self.r_matrix = kalman_coef['r']
        self.orientation = orientation
        self.single_femto_filter = single_femto_filter
        self.sequence_under_diagonal = flatten([[[i, j] for i in range(j)] for j in range(self.f.n)])
        if not self.orientation:
            self.r_orf_estimation = np.array([f.rv_orf_calc[i][0:6] for i in range(f.n)])
            self.phi_ = np.array([[1, 0, 0, dt, 0, 0],
                                  [0, 1, 0, 0, dt, 0],
                                  [0, 0, 1, 0, 0, dt],
                                  [0, 0, 0, 1, 0, -2 * w_orb * dt],
                                  [0, - w_orb ** 2 * dt, 0, 0, 1, 0],
                                  [0, 0, 3 * w_orb ** 2 * dt,  2 * w_orb * dt, 0, 1]])
            self.d_ = np.vstack([np.zeros((3, 3)), np.eye(3)])
            self.p_ = [np.diag([self.p_coef[0]]*3 + [self.p_coef[1]]*3) for _ in range(f.n)]
        else:  # 'rvq'
            self.r_orf_estimation = np.array(f.rv_orf_calc)
            self.phi_ = np.array([[1, 0, 0, dt, 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, dt, 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, dt, 0, 0, 0],
                                  [0, 0, 0, 1, 0, -2 * w_orb * dt, 0, 0, 0],
                                  [0, -w_orb ** 2 * dt, 0, 0, 1, 0, 0, 0, 0],
                                  [0, 0, 3 * w_orb ** 2 * dt,  2 * w_orb * dt, 0, 1, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 1, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 1]])
            self.d_ = np.vstack([np.zeros((3, 3)), np.eye(3), np.eye(3)])
            self.p_ = [np.diag([self.p_coef[0]]*3 + [self.p_coef[1]]*3 + [self.p_coef[2]]*3) for _ in range(f.n)]
        self.q_ = np.eye(3) * self.q_coef

        t = 9 if self.orientation else 6
        self.sigmas = [[] for _ in range(t * self.f.n)]
        self.real_sigmas = [[] for _ in range(t * self.f.n)]
        if not self.single_femto_filter:  # Расширешние на учёт несколько аппаратов а фильтре
            self.phi_ = np.bmat([[np.zeros([t, t])] * i + [self.phi_] +
                                 [np.zeros([t, t])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.q_ = np.bmat([[np.zeros([3, 3])] * i + [self.q_] +
                               [np.zeros([3, 3])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.p_ = np.bmat([[np.zeros([t, t])] * i + [self.p_[0]] +
                               [np.zeros([t, t])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.d_ = np.bmat([[np.zeros([t, 3])] * i + [self.d_] +
                               [np.zeros([t, 3])] * (self.f.n - i - 1) for i in range(self.f.n)])
            # self.d_ = np.bmat([[self.d_] for i in range(self.f.n)])

    def calc(self, i: int) -> None:
        t = 9 if self.orientation else 6
        z_ = self.c.calc_dist[0][i][-1]
        r_ = self.r_orf_estimation[i]
        if self.p.is_aero:
            rv_m = self.p.aero_integrate(self.f, i=i, r=r_[0:3], v=r_[3:6])
            r_m = np.append(rv_m[0], rv_m[1]) if t == 6 else np.append(np.append(rv_m[0], rv_m[1]), np.zeros(3))
        else:
            r_m = self.phi_ @ r_

        if self.c.gain_mode == self.c.gain_modes[4]:
            pass
            '''if self.orientation:
                tmp1 = self.c.get_gain(quart2dcm(vec2quat(np.zeros(3))) @ np.array(self.c.r_orf[i] - r_m[0:3]))
                tmp2 = self.f.get_gain(quart2dcm(vec2quat(r_m[6:9])) @ np.array(self.c.r_orf[i] - r_m[0:3]), mode3=True)
                signal_rate = [tmp1[ii] * tmp2 for ii in range(2)]
                z_model = np.array([np.linalg.norm(r_m[0:3] - self.c.r_orf[0]) / np.sqrt(signal_rate[i])
                                    for i in range(2)])
            else:
                z_model = np.array([1., 1.])
            z_ = self.c.calc_dist_[0][i]

            h_ = np.array([[(r_m[j] - self.c.r_orf[0][j]) / z_model[ii] for j in range(3)] + (t - 3) * [0.]
                           for ii in range(2)])
            q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt
            p_m = self.phi_ @ self.p_[i] @ self.phi_.T + q_tilda
            k_ = p_m @ h_.T @ np.linalg.pinv(h_ @ p_m @ h_.T + np.eye(2) * self.r_matrix)
            self.r_orf_estimation[i] = r_m + k_ @ (z_ - z_model)
            self.p_[i] = (np.eye(t) - k_ @ h_) @ p_m'''
        else:
            signal_rate = self.c.get_gain(np.array(self.c.r_orf[i] - r_m[0:3])) * \
                self.f.get_gain(quart2dcm(vec2quat(r_m[6:9])) @ np.array(self.c.r_orf[i] - r_m[0:3])) \
                if self.orientation else 1
            z_model = np.linalg.norm(r_m[0:3] - self.c.r_orf[0])
            z_ *= np.sqrt(signal_rate)
            self.f.z_difference[i] += [abs(z_model - z_)]

            h_ = np.array([(r_m[j] - self.c.r_orf[0][j]) / z_model for j in range(3)] + (t - 3) * [0.])
            q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt
            p_m = self.phi_ @ self.p_[i] @ self.phi_.T + q_tilda
            k_ = p_m @ h_.T / (h_ @ p_m @ h_.T + self.r_matrix)
            self.r_orf_estimation[i] = r_m + k_ * (z_ - z_model)
            self.p_[i] = (np.eye(t) - np.outer(k_, h_)) @ p_m

    def calc_all(self) -> None:
        t = 9 if self.orientation else 6
        ll = int(self.f.n*(self.f.n + 1)/2)  # Рёбра полного графа N фемтосатов + 1 кубсат
        z_ = [self.c.calc_dist[0][i][-1] for i in range(self.f.n)] + flatten([
            [self.f.calc_dist[j][i][-1] for i in range(j)] for j in range(self.f.n)])
        if self.p.is_aero:
            r_ = self.r_orf_estimation
            rv_m = [self.p.aero_integrate(self.f, i=i, r=r_[i][0:3], v=r_[i][3:6]) for i in range(self.f.n)]
            r_m = np.array(flatten([np.append(rv_m[i][0], rv_m[i][1]) if t == 6 else
                                    np.append(np.append(rv_m[i][0], rv_m[i][1]), np.zeros(3))
                                    for i in range(self.f.n)]))
        else:
            r_ = np.array(flatten(self.r_orf_estimation))
            r_m = self.phi_ @ r_
            r_m = np.array(self.phi_ @ r_)[0]  # 6t_6t @ 6t -> 6t

        if self.orientation:  # ВЕРНУТЬСЯ СЮДА
            signal_rate = np.array([self.f.get_gain(quart2dcm(vec2quat(r_m[9*i+6:9*i+9])) @
                                                    (r_m[9*i+0:9*i+3] - self.c.r_orf[0])) *
                                    self.c.get_gain(r_m[9*i+0:9*i+3] - self.c.r_orf[0])
                                    for i in range(self.f.n)] +
                                   flatten([[self.f.get_gain(quart2dcm(vec2quat(r_m[9*i+6:9*i+9]))
                                                             @ (r_m[9*i+0:9*i+3] - r_m[9*j+0:9*j+3])) *
                                             self.f.get_gain(quart2dcm(vec2quat(r_m[9*j+6:9*j+9]))
                                                             @ (r_m[9*i+0:9*i+3] - r_m[9*j+0:9*j+3]))
                                             for i in range(j)] for j in range(self.f.n)]))
        else:
            signal_rate = np.array([1] * int(self.f.n * (self.f.n + 1) // 2))
        z_model = np.array([np.linalg.norm(r_m[0+t*i:3+t*i] - self.c.r_orf[0]) for i in range(self.f.n)] + flatten([
            [np.linalg.norm(r_m[0+t*j:3+t*j] - r_m[0+t*i:3+t*i]) for i in range(j)] for j in range(self.f.n)]))
        z_ *= np.sqrt(signal_rate)
        for i in range(self.f.n):
            self.f.z_difference[i] += [abs(z_model[i] - z_[i])]

        def local_h_func(i: int, j: int) -> list:
            """Функция для построения матрицы H:
            i - строка из n + n(n-1)/2
            j - столбец из n (6 строк выдаёт как блок в return)"""
            r_k = r_m[t*j:t*j+3]
            if (i < self.f.n and i == j) or (i >= self.f.n and j in self.sequence_under_diagonal[i - self.f.n]):
                return [(r_k[j] - self.c.r_orf[0][j]) / z_model[i] for j in range(3)] + (t - 3) * [0.]
            else:
                return t * [0.]

        h_ = np.vstack([flatten([local_h_func(i, j) for i in range(self.f.n)]) for j in range(ll)])
        q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt  # nt_nt
        p_m = self.phi_ @ self.p_ @ self.phi_.T + q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        r = np.eye(ll) * self.r_matrix  # n + n(n-1)/2
        k_ = p_m @ h_.T @ np.linalg.inv(h_ @ p_m @ h_.T + r)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.p_ = (np.eye(t * self.f.n) - k_ @ h_) @ p_m
        r_orf_estimation = np.matrix(r_m) + k_ @ (z_ - z_model)
        for i in range(self.f.n):
            self.r_orf_estimation[i] = np.array(r_orf_estimation)[0][0+i*t:t+i*t]

class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat, s: Cosmetic, kalman_coef: dict, method_navigation: str, dt: float = 1.,
                 is_aero: bool = True, is_complex_aero: bool = False, is_hkw: bool = True, h_orb: float = 600e3):
        self.dt = dt
        self.h_orb = h_orb
        self.show_rate = 1

        # Неизменные параметры
        self.t = 0.
        self.iter = 0
        self.mu = 5.972e24 * 6.67408e-11  # гравитационный параметр
        self.j_2 = 1.082 * 1e-3
        self.r_earth = 6371e3
        self.r_orb = self.r_earth + self.h_orb
        self.w_orb = np.sqrt(self.mu / self.r_orb ** 3)
        self.v_orb = np.sqrt(self.mu / self.r_orb)
        self.s = s
        self.c = c
        self.f = f
        self.r_matrix = kalman_coef['r']
        self.time_begin = datetime.now()

        # Параметры типа bool
        self.is_aero = is_aero
        self.is_complex_aero = is_complex_aero
        self.is_hkw = is_hkw

        # Параметры фильтров
        self.navigation = method_navigation.split()
        self.k = KalmanFilter(f=f, c=c, s=s, p=self, dt=self.dt, w_orb=self.w_orb, kalman_coef=kalman_coef,
                              orientation='rvq' in self.navigation, single_femto_filter='all' not in self.navigation)

        # Подгон чтобы не разлеталися (C_1=0 у всех)
        self.c.c_hkw = [get_c_hkw(self.c.r_orf[i], self.c.v_orf[i], self.w_orb) for i in range(self.c.n)]
        self.f.c_hkw = [get_c_hkw(self.f.r_orf[i], self.f.v_orf[i], self.w_orb) for i in range(self.f.n)]

    # Интегрирование движения
    def aero_integrate(self, obj, i: int, r=None, v=None):
        S = quart2dcm(obj.q[i])
        cos_alpha = clip((np.trace(S) - 1) / 2, -1, 1)
        # alpha = 180 / np.pi * np.arccos(cos_alpha)
        rho = self.get_atm_params(obj.r_orf[i][2])[0]
        if obj is self.c:
            c_resist = 1.05
            square = obj.size[0] * obj.size[1]
        else:
            c_resist = 1.17
            square = obj.size[0] * obj.size[1] * abs(cos_alpha)
        r = obj.r_orf[i] if r is None else r
        v = obj.v_orf[i] if v is None else v
        return self.rk4_translate(r, v,
                                  self.get_full_acceleration(c_resist=c_resist, rho=rho, square=square, m=obj.mass,
                                                             r=r, v=v))

    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt

        # Утверждение параметров
        if self.iter == 1:
            tmp1 = f", сферичность={int(self.c.ellipsoidal_signal*100)}%" if self.c.gain_mode == self.c.gain_modes[1] \
                else ""
            tmp2 = f", сферичность={int(self.f.ellipsoidal_signal*100)}%" if self.f.gain_mode == self.f.gain_modes[1] \
                else ""
            tmp = ", ориентации\n" if 'rvq' in self.navigation else "\n"
            self.s.my_print(f"Диаграмма антенн кубсата: {self.c.gain_mode}{tmp1}\n"
                            f"Диаграмма антенн фемтосатов: {self.f.gain_mode}{tmp2}\n"
                            f"Учёт аэродинамики: {self.is_aero}\n"
                            f"Применяется фильтр Калмана для поправки: положений, скоростей{tmp}" 
                            f"Фильтр Калмана основан на: "
                            f"{'одном чипсате' if self.k.single_femto_filter else 'всех чипсатах'}", color='c')

        # Вращательное движение
        # -----его нет

        # Поступательное движение
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                if self.is_aero:
                    obj.r_orf[i], obj.v_orf[i] = self.aero_integrate(obj, i)
                else:
                    obj.r_orf[i] = r_hkw(obj.c_hkw[i], self.w_orb, self.t)
                    obj.v_orf[i] = v_hkw(obj.c_hkw[i], self.w_orb, self.t)

                if self.iter % self.show_rate == 0:
                    obj.line[i] += [obj.r_orf[i][0], obj.r_orf[i][1], obj.r_orf[i][2]]

        # Расчёт связи
        for i_c in range(self.c.n):
            # S_c = quart2dcm(self.c.q[i_c])
            for i_f in range(self.f.n):
                r_tmp = self.f.r_orf[i_f] - self.c.r_orf[i_c]
                self.c.real_dist[i_c][i_f] += [np.linalg.norm(r_tmp)]
                '''
                S_f = quart2dcm(self.f.q[i_f])
                signal = self.c.get_gain(r=S_c @ r_tmp, mode3=True) * self.f.get_gain(r=S_f @ r_tmp, mode3=True) \
                #     * self.f.power_signal_full / self.f.length_signal_full**2 / np.linalg.norm(r_tmp)**2
                # signal *= 1 + np.random.normal(0, np.sqrt(self.r_matrix))
                signal += np.random.normal(0, np.sqrt(self.r_matrix)) * \
                    self.f.power_signal_full / self.f.length_signal_full**2 / np.linalg.norm(r_tmp)**2
                # self.c.signal_power[i_c][i_f] += [signal]'''
                S_tmp = quart2dcm(np.array(self.k.f.q[i_f]))
                calc_dist = np.linalg.norm(r_tmp) / np.sqrt(self.f.get_gain(r=S_tmp @ r_tmp, mode3=True) *
                                                            self.c.get_gain(r=r_tmp, mode3=True)) + \
                    np.random.normal(0, np.sqrt(self.r_matrix))
                """if self.c.gain_mode == self.c.gain_modes[4]:
                    pass
                    tmp = self.f.get_gain(r=S_tmp @ r_tmp)
                    calc_dist_ = np.array([np.sqrt(self.f.power_signal_full/self.f.length_signal_full**2 / signal
                                                   * tmp[ii]) for ii in range(2)])
                calc_dist *= 1 + np.random.normal(0, np.sqrt(self.r_matrix))
                if self.c.gain_mode == self.c.gain_modes[4]:
                    self.c.calc_dist_[i_c][i_f] = calc_dist_"""
                self.c.calc_dist[i_c][i_f] += [calc_dist]
                self.c.kalm_dist[i_c][i_f] += [np.linalg.norm(self.f.r_orf[i_f] -
                                                              self.k.r_orf_estimation[i_f][0:3])]
                if self.iter % self.show_rate == 0:
                    self.f.line_kalman[i_f] += [self.k.r_orf_estimation[i_f][0],
                                                self.k.r_orf_estimation[i_f][1],
                                                self.k.r_orf_estimation[i_f][2]]
                    self.f.line_difference[i_f] += \
                        [np.linalg.norm(self.k.r_orf_estimation[i_f][0:3] - np.array(self.f.r_orf[i_f]))]
        for i_f1 in range(self.f.n):
            for i_f2 in range(self.f.n):
                if i_f1 != i_f2:
                    r_tmp = self.f.r_orf[i_f1] - self.f.r_orf[i_f2]
                    self.f.real_dist[i_f1][i_f2] += [np.linalg.norm(r_tmp)]

                    q1 = self.k.f.q[i_f1]
                    q2 = self.k.f.q[i_f1]
                    calc_dist = np.linalg.norm(r_tmp) / np.sqrt(self.f.get_gain(r=quart2dcm(q1) @ r_tmp, mode3=True) *
                                                                self.f.get_gain(r=quart2dcm(q2) @ r_tmp, mode3=True)) +\
                        np.random.normal(0, np.sqrt(self.r_matrix))
                    self.f.calc_dist[i_f1][i_f2] += [calc_dist]
                else:
                    self.f.calc_dist[i_f1][i_f2] += [0.]

        # Оценка положения фемтоспутников
        t = 9 if self.k.orientation else 6
        if self.k.single_femto_filter:
            for i in range(self.f.n):
                self.k.calc(i)
                if self.c.gain_mode != self.c.gain_modes[4]:  # 3
                    for j_n in range(self.f.n):
                        for j_t in range(t):
                            tmp = np.append(np.append(self.f.r_orf[j_n], self.f.v_orf[j_n]), list(self.f.q[j_n][1:4]))
                            self.k.sigmas[j_n * t + j_t] += [np.sqrt(self.k.p_[j_n][j_t][j_t]) * tmp[j_t]]
        else:
            self.k.calc_all()
            '''if False:  # self.c.gain_mode != self.c.gain_modes[4]:  # грёбаных
                for j in range(t * self.f.n):
                    self.k.sigmas[j] += [np.sqrt(self.k.p_[j][j])]
        if False:  # self.c.gain_mode != self.c.gain_modes[4]:  # раза
            for j_n in range(self.f.n):
                for j_t in range(t):
                    tmp = np.append(np.append(self.f.r_orf[j_n], self.f.v_orf[j_n]), list(self.f.q[j_n][1:4]))
                    self.k.real_sigmas[j_n * t + j_t] += [np.abs(tmp - np.array(self.k.r_orf_estimation[j_n][j_t]))]'''

    def integrate(self, t: float):
        self.s.my_print(f"Оборотов вокруг Земли: {round(10 * t / (3600 * 1.5)) / 10}\n"
                        f"Дней: {round(100 * t / (3600 * 24)) / 100}", color='b')
        n = int(t//self.dt)
        flag = 0.
        for i in range(n):
            if i / n > (flag + 0.1):
                flag += 0.1
                per = int(10 * i / n)
                self.s.my_print(f"{10 * per}% [{'#' * per + ' ' * (10 - per)}]"
                                + real_workload_time(n=per, n_total=10, time_begin=self.time_begin,
                                                     time_now=datetime.now()), color='m')
                # print(f"P:{[np.linalg.det(self.k.p_[ii]) for ii in range(self.f.n)]}")
            self.time_step()

    # Функции возврата
    def get_full_acceleration(self, c_resist: float, rho: float, square: float, r: Union[float, np.ndarray],
                              v: Union[float, np.ndarray], m: float) -> np.ndarray:
        """Возвращает вектор силы в ОСК, принимает параметры в ОСК\n
        square - площадь S для аэродинамики
        c_resist - Cf для аэродинаммики"""
        # force = - r * self.mu / np.linalg.norm(r)
        force = self.get_w_orb_acceleration(r, v)
        v_real = v + np.array([self.v_orb, 0, 0])
        if self.is_aero:
            force -= v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / m
        return force

    def get_w_orb_acceleration(self, r, v):
        return np.array([-2 * self.w_orb * v[2],
                         -self.w_orb ** 2 * r[1],
                         2 * self.w_orb * v[0] + 3 * self.w_orb ** 2 * r[2]])

    def get_atm_params(self, h: float):
        """https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html"""
        t = -131.21 + 0.00299 * (h + self.h_orb)
        p = 2.488 * ((t + 273.1) / 216.6) ** -11.388
        rho = p / (0.2869 * (t + 273.1))
        return rho, t, p

    # Рунге-Кутты 4 порядка
    def rk4_translate(self, r, v, a):
        def rv_right_part(rv1, a1):
            return np.array([rv1[3], rv1[4], rv1[5], a1[0], a1[1], a1[2]])
        rv = np.append(r, v)
        k1 = rv_right_part(rv, a)
        k2 = rv_right_part(rv + k1 * self.dt / 2, a)
        k3 = rv_right_part(rv + k2 * self.dt / 2, a)
        k4 = rv_right_part(rv + k3 * self.dt, a)
        rv = self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return rv[0:3] + r, rv[3:6] + v

class Objects:
    def __init__(self, method_navigation: str, kalman_coef: dict, dt: float = 1, n_c: int = 1, n_f: int = 5,
                 model_c: str = '1U', if_any_print: bool = True, start_navigation: str = 'perfect',
                 start_navigation_tolerance: float = 0.9):
        """Класс содержит информацию о n_c кубсатах и n_f фемтосатах. \n
        Размер n_c кубсатов определяется моделью model_c"""
        # Проверки на вшивость
        models = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        if model_c not in models:
            raise ValueError("Внимание! Модель кубсата должна быть одна из: 1U/1.5U/2U/3U/6U/12U")

        self.n_a = 1
        self.n_f = n_f
        self.model_c = model_c

        # Классы
        self.h_orb = 600e3
        self.mu = 5.972e24 * 6.67408e-11  # гравитационный параметр
        self.j_2 = 1.082 * 1e-3
        self.r_earth = 6371e3
        self.r_orb = self.r_earth + self.h_orb
        self.w_orb = np.sqrt(self.mu / self.r_orb ** 3)
        self.s = Cosmetic(if_any_print=if_any_print)
        self.c = CubeSat(n=n_c, n_f=n_f, model=model_c, w_orb=self.w_orb)
        self.f = FemtoSat(n=n_f, start_navigation_tolerance=start_navigation_tolerance,
                          start_navigation=start_navigation, w_orb=self.w_orb)
        self.p = PhysicModel(c=self.c, f=self.f, s=self.s, dt=dt, method_navigation=method_navigation,
                             kalman_coef=kalman_coef, h_orb=self.h_orb)
