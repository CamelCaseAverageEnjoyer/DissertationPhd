from colorama import init, Fore, Style, Back
from typing import Union, Any, Callable
from random import uniform
import scipy
import random
import numpy as np
import kiam_astro

from tiny_functions import *


class FemtoSat:
    def __init__(self, n: int = 10, r_spread: float = 1e2, v_spread: float = 1e-1):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        # Предопределённые параметры
        self.gain_modes = ['isotropic', 'ellipsoid']

        # Общие параметры
        self.n = n
        self.size = [0.03, 0.03]
        self.mass = 0.1
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001

        self.ellipsoidal_signal = 0.9
        self.gain_mode = self.gain_modes[1]

        # Индивидуальные параметры
        self.r_orf_calc = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.r_orf = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.v_orf = [np.array([uniform(-v_spread, v_spread) for _ in range(3)]) for _ in range(self.n)]
        self.c_hkw = [[] for i in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]
        self.signal_power = [[[] for _ in range(self.n)] for _ in range(self.n)]
        self.real_dist = [[[] for _ in range(self.n)] for _ in range(self.n)]
        self.calc_dist = [[[] for _ in range(self.n)] for _ in range(self.n)]

    def get_gain(self, r: Union[float, np.ndarray]):
        r1 = r / np.linalg.norm(r)
        if self.gain_mode == self.gain_modes[1]:
            return np.linalg.norm([r1[0] * 1, r1[1] * self.ellipsoidal_signal, r1[2] * self.ellipsoidal_signal])
        return 1

class CubeSat:
    def __init__(self, n_f: int, n: int = 1, model: str = '1U', r_spread: float = 100, v_spread: float = 0.1):
        """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.\n
        Все величны представлены в СИ."""
        # Предопределённые параметры
        self.gain_modes = ['isotropic', 'ellipsoid']
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
        self.c_hkw = [[] for i in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]
        self.signal_power = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.real_dist = [[[] for _ in range(n_f)] for _ in range(self.n)]
        self.calc_dist = [[[] for _ in range(n_f)] for _ in range(self.n)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

    def get_gain(self, r: Union[float, np.ndarray]):
        r1 = r / np.linalg.norm(r)
        if self.gain_mode == self.gain_modes[1]:
            return np.linalg.norm([r1[0] * 1, r1[1] * self.ellipsoidal_signal, r1[2] * self.ellipsoidal_signal])
        return 1

class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat, dt: float = 1., is_aero: bool = False, is_complex_aero: bool = False,
                 is_hkw: bool = True):
        self.dt = dt
        self.h_orb = 600e3
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
        self.c = c
        self.f = f

        # Параметры филтров
        self.r_orf_estimation = np.array(self.f.r_orf_calc)
        self.phi_ = np.array([[1, 0, 0, self.dt, 0, 0],
                              [0, 1, 0, 0, self.dt, 0],
                              [0, 0, 1, 0, 0, self.dt],
                              [0, 0, 0, 1, 0, -2 * self.w_orb * self.dt],
                              [0, -self.w_orb ** 2 * self.dt, 0, 0, 1, 0],
                              [0, 0, 3 * self.w_orb ** 2 * self.dt,  2 * self.w_orb * self.dt, 0, 1]])
        self.q_ = np.eye(6) * 1e-3
        self.p_ = np.eye(6) * 1e-3
        self.r_ = 1e-3

        # Параметры типа bool
        self.is_aero = is_aero
        self.is_complex_aero = is_complex_aero
        self.is_hkw = is_hkw

        # Подгон чтобы не разлеталися
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                obj.q[i] /= np.linalg.norm(obj.q[i])
                obj.v_orf[i][0] = - 2 * obj.r_orf[i][2] * self.w_orb
        self.c.c_hkw = [get_c_hkw(self.c.r_orf[i], self.c.v_orf[i], self.w_orb) for i in range(self.c.n)]  # Рудимент
        self.f.c_hkw = [get_c_hkw(self.f.r_orf[i], self.f.v_orf[i], self.w_orb) for i in range(self.f.n)]  # Рудимент

    # Интегрирование движения
    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt

        # Утверждение параметров
        if self.iter == 1:
            tmp1 = f", сферичность={int(self.c.ellipsoidal_signal*100)}%" if self.c.gain_mode == self.c.gain_modes[1] \
                else ""
            tmp2 = f", сферичность={int(self.f.ellipsoidal_signal*100)}%" if self.f.gain_mode == self.f.gain_modes[1] \
                else ""
            print(f"Диаграмма антенн кубсата: {self.c.gain_mode}{tmp1}")
            print(f"Диаграмма антенн фемтосатов: {self.f.gain_mode}{tmp2}")

        # Вращательное движение
        # -----его нет

        # Поступательное движение
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                if self.is_aero:
                    S = quart2dcm(obj.q[i])
                    cos_alpha = clip((np.trace(S) - 1) / 2, -1, 1)
                    # alpha = 180 / np.pi * np.arccos(cos_alpha)
                    rho = self.get_atm_params(self.h_orb)[0]
                    if obj is self.c:
                        c_resist = 1.05
                        square = obj.size[0] * obj.size[1]
                    else:
                        c_resist = 1.17
                        square = obj.size[0] * obj.size[1] * abs(cos_alpha)
                    obj.r_orf[i], obj.v_orf[i] = self.rk4_translate(obj.r_orf[i], obj.v_orf[i],
                                                                    self.get_full_acceleration(c_resist=c_resist,
                                                                                               rho=rho,
                                                                                               square=square,
                                                                                               m=obj.mass,
                                                                                               r=obj.r_orf[i],
                                                                                               v=obj.v_orf[i]))
                else:
                    obj.r_orf[i] = r_hkw(obj.c_hkw[i], self.w_orb, self.t)
                    obj.v_orf[i] = v_hkw(obj.c_hkw[i], self.w_orb, self.t)

                if self.iter % self.show_rate == 0:
                    obj.line[i] += [obj.r_orf[i][0], obj.r_orf[i][1], obj.r_orf[i][2]]

        # Расчёт связи
        if self.iter % self.show_rate == 0:
            for i_c in range(self.c.n):
                S_c = quart2dcm(self.c.q[i_c])
                for i_f in range(self.f.n):
                    S_f = quart2dcm(self.f.q[i_f])
                    r_tmp = self.f.r_orf[i_f] - self.c.r_orf[i_c]
                    signal = self.c.get_gain(r=S_c @ r_tmp) * self.f.get_gain(r=S_f @ r_tmp) \
                        * self.f.power_signal_full / self.f.length_signal_full**2 / np.linalg.norm(r_tmp)**2
                    self.c.real_dist[i_c][i_f] += [np.linalg.norm(r_tmp)]
                    self.c.signal_power[i_c][i_f] += [signal]
                    self.c.calc_dist[i_c][i_f] += [np.sqrt(self.f.power_signal_full/self.f.length_signal_full**2
                                                           / signal)]
            for i_f1 in range(self.f.n):
                S_f1 = quart2dcm(self.f.q[i_f1])
                for i_f2 in range(self.f.n):
                    if i_f1 != i_f2:
                        S_f2 = quart2dcm(self.f.q[i_f2])
                        r_tmp = self.f.r_orf[i_f1] - self.f.r_orf[i_f2]
                        signal = self.f.get_gain(r=S_f1 @ r_tmp) * self.f.get_gain(r=S_f2 @ r_tmp) * \
                            self.f.power_signal_full / self.f.length_signal_full**2 / np.linalg.norm(r_tmp)**2
                        self.f.real_dist[i_f1][i_f2] += [np.linalg.norm(r_tmp)]
                        self.f.signal_power[i_f1][i_f2] += [signal]
                        self.f.calc_dist[i_f1][i_f2] += [np.sqrt(self.f.power_signal_full/self.f.length_signal_full**2
                                                                 / signal)]

        # Оценка положения фемтоспутников
        for i in range(self.f.n):
            d_ = self.c.calc_dist[0][i][len(self.c.calc_dist[0][i]) - 1]
            r_ = self.r_orf_estimation[i]
            h_ = np.array([r_[0]/d_, r_[1]/d_, r_[2]/d_])
            self.q_ = self.phi_ @ self.q_ @ self.phi_.T * self.dt
            print(f"q_={np.linalg.det(self.q_)}")

            r_m = self.phi_ @ r_
            p_m = self.phi_ @ self.p_ @ self.phi_.T + self.q_

            k_ = p_m @ h_.T / (h_ @ p_m @ h_.T + self.r_)
            self.r_orf_estimation = r_m + k_ @ (1 - h_ * r_m)
            self.p_ = (np.eye(6) - k_ * h_) @ p_m
        '''a_m_f = np.zeros((self.f.n, self.f.n))
        a_p_f = np.zeros((self.f.n, self.f.n))
        a_m_c = np.zeros((self.f.n, self.c.n))
        a_p_c = np.zeros((self.f.n, self.c.n))

        def local_func(x1):
            anw = 0
            for i__c in range(self.c.n):
                for i__f in range(self.f.n):
                    anw += a_m_c[i__f][i__c] + a_p_c[i__f][i__c]
            for i_1f in range(self.c.n):
                for i_2f in range(self.f.n):
                    anw += a_m_f[i_1f][i_2f] + a_p_f[i_1f][i_2f]
            return anw

        def local_constraints(x1):
            anw = []
            x = x1.reshape((self.f.n, 3))
            z = np.vstack([np.hstack((np.eye(3), x.T)), np.hstack((x, x @ x.T))])

            for i__c in range(self.c.n):
                for i__f in range(self.f.n):
                    e = np.hstack([[int(j == i__f) for j in range(self.f.n)], self.c.r_orf[i__c]])
                    anw += [-self.c.calc_dist[i__c][i__f][len(self.c.calc_dist[i__c][i__f]) - 1] ** 2 +
                            e.T @ z @ e - a_p_c[i__f][i__c] + a_m_c[i__f][i__c]]
            for i_1f in range(self.c.n):
                for i_2f in range(self.f.n):
                    if i_1f != i_2f:
                        e = np.hstack([[int(j == i_1f) - int(j == i_2f) for j in range(self.f.n)], np.zeros(3)])
                        anw += [-self.f.calc_dist[i_1f][i_2f][len(self.f.calc_dist[i_1f][i_2f]) - 1] ** 2 +
                                e.T @ z @ e - a_p_f[i_1f][i_2f] + a_m_f[i_1f][i_2f]]
            return anw

        tol = 1
        opt = {'verbose': 3, 'gtol': 1e-20, 'xtol': 1e-20}
        x0 = flatten(self.r_orf_estimation)
        res = scipy.optimize.minimize(local_func, x0, tol=tol, method='trust-constr', options=opt,
                                      constraints={'type': 'eq', 'fun': local_constraints})
        self.r_orf_estimation -= 0'''

    def integrate(self, t: float):
        for _ in range(int(t//self.dt)):
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
            # force += - v * np.linalg.norm(v) * c_resist * rho * square / 2
            force -= v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / m
        return force

    def get_w_orb_acceleration(self, r, v):
        return np.array([-2 * self.w_orb * v[2],
                         -self.w_orb ** 2 * r[1],
                         2 * self.w_orb * v[0] + 3 * self.w_orb ** 2 * r[2]])

    def get_atm_params(self, h: float):
        """https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html"""
        t = -131.21 + 0.00299 * h
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
    def __init__(self, dt: float = 1, n_c: int = 1, n_f: int = 5, model_c: str = '1U', if_any_print: bool = True):
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
        self.c = CubeSat(n=n_c, n_f=n_f, model=model_c)
        self.f = FemtoSat(n=n_f)
        self.p = PhysicModel(c=self.c, f=self.f, dt=dt)

        # Косметика
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
