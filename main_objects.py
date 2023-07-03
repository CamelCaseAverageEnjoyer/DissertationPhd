from colorama import init, Fore, Style, Back
from typing import Union, Any, Callable
from random import uniform
import random
import numpy as np
import kiam_astro

from tiny_functions import *


class FemtoSat:
    def __init__(self, n: int = 10, r_spread: float = 1e2, v_spread: float = 1e-1):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""

        # Общие параметры
        self.n = n
        self.size = [0.03, 0.03]
        self.mass = 0.1
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001

        # Индивидуальные параметры
        self.r_orf = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.v_orf = [np.array([uniform(-v_spread, v_spread) for _ in range(3)]) for _ in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]

    def get_gain(self, r: Union[float, np.ndarray], mode: str = 'isotropic'):
        gain_modes = ['ellipsoid']
        if mode not in gain_modes:
            mode = 'isotropic'
        # print(f"Диаграмма антенн фемтосатов: {mode}")
        r1 = r / np.linalg.norm(r)
        if mode == gain_modes[0]:
            return np.linalg.norm([r1[0] * 1, r1[1] * 0.5, r1[2] * 0.5])
        return 1

class CubeSat:
    def __init__(self, n: int = 1, model: str = '1U', r_spread: float = 100, v_spread: float = 0.1):
        """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.\n
        Все величны представлены в СИ."""
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

        # Индивидуальные параметры
        self.r_orf = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.v_orf = [np.array([uniform(-v_spread, v_spread) for _ in range(3)]) for _ in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]
        self.signal_power = [[[] for _ in range(self.n)] for _ in range(self.n)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

    def get_gain(self, r: Union[float, np.ndarray], mode: str = 'isotropic'):
        gain_modes = ['ellipsoid']
        if mode not in gain_modes:
            mode = 'isotropic'
        # print(f"Диаграмма антенн кубсатов: {mode}")
        r1 = r / np.linalg.norm(r)
        if mode == gain_modes[0]:
            return np.linalg.norm([r1[0] * 1, r1[1] * 0.5, r1[2] * 0.5])
        return 1

class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat, dt: float = 1., is_aero: bool = True, is_complex_aero: bool = False,
                 is_hkw: bool = True, show_rate: int = 1):
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
        self.C_c = [get_c_hkw(self.c.r_orf[i], self.c.v_orf[i], self.w_orb) for i in range(self.c.n)]  # Рудимент
        self.C_f = [get_c_hkw(self.f.r_orf[i], self.f.v_orf[i], self.w_orb) for i in range(self.f.n)]  # Рудимент
        self.c.signal_power = [[[] for _ in range(self.f.n)] for _ in range(self.c.n)]

        # Параметры типа bool
        self.is_aero = is_aero
        self.is_complex_aero = is_complex_aero
        self.is_hkw = is_hkw

        # Подгон
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                obj.q[i] /= np.linalg.norm(obj.q[i])
                obj.v_orf[i][0] = - 2 * obj.r_orf[i][2] * self.w_orb

    # Интегрирование движения
    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt

        # Вращательное движение
        # -----похуй
        # Поступательное движение
        for obj in [self.c, self.f]:
            for i in range(obj.n):
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
                                                                self.get_full_acceleration(c_resist=c_resist, rho=rho,
                                                                                           square=square, m=obj.mass,
                                                                                           r=obj.r_orf[i],
                                                                                           v=obj.v_orf[i]))
                if self.iter % self.show_rate == 0:
                    obj.line[i] += [obj.r_orf[i][0], obj.r_orf[i][1], obj.r_orf[i][2]]

        # Связь
        if self.iter % self.show_rate == 0:
            for i_c in range(self.c.n):
                S_c = quart2dcm(self.c.q[i_c])
                for i_f in range(self.f.n):
                    S_f = quart2dcm(self.f.q[i_f])
                    r_tmp = self.f.r_orf[i_c] - self.c.r_orf[i_c]
                    self.c.signal_power[i_c][i_f] += [self.c.get_gain(r=S_c @ r_tmp, mode='isotropic') *
                                                      self.f.get_gain(r=S_f @ r_tmp, mode='ellipsoid') *
                                                      self.f.power_signal_full / self.f.length_signal_full**2 /
                                                      np.linalg.norm(r_tmp)**2]

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
    def rv_right_part(self, rv, a):
        return np.array([rv[3], rv[4], rv[5], a[0], a[1], a[2]])

    def rk4_translate(self, r, v, a):
        rv = np.append(r, v)
        k1 = self.rv_right_part(rv, a)
        k2 = self.rv_right_part(rv + k1 * self.dt / 2, a)
        k3 = self.rv_right_part(rv + k2 * self.dt / 2, a)
        k4 = self.rv_right_part(rv + k3 * self.dt, a)
        rv = self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return rv[0:3] + r, rv[3:6] + v

class Objects:
    def __init__(self, dt: float = 1, n_c: int = 1, n_f: int = 5, model_c: str = '1U', if_any_print: bool = True,
                 show_rate: int = 1):
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
        self.c = CubeSat(n=n_c, model=model_c)
        self.f = FemtoSat(n=n_f)
        self.p = PhysicModel(c=self.c, f=self.f, dt=dt, show_rate=show_rate)

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
