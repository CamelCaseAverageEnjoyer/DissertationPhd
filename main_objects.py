from colorama import init, Fore, Style, Back
from typing import Union, Any, Callable
from random import uniform
import random
import numpy as np
import kiam_astro as ka


class PhysicModel:
    def __init__(self, is_aero: bool = True, is_complex_aero: bool = False, is_j_2: bool = True):
        self.mu = 5.972e24 * 6.67408e-11  # гравитационный параметр
        self.is_aero = is_aero
        self.is_complex_aero = is_complex_aero
        self.j_2 = 1.082 * 1e-3

    def get_force(self, c_resist: float, rho: float, square: float, r: Union[float, np.ndarray],
                  v: Union[float, np.ndarray]) -> np.ndarray:
        """Возвращает вектор силы в ИСК, принимает параметры в ИСК"""
        force = - r * self.mu / np.linalg.norm(r)
        if self.is_aero:
            force += - v * np.linalg.norm(v) * c_resist * rho * square / 2
        return force

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
        self.size = sizes[self.model_number]

        # Индивидуальные параметры
        self.r = [np.array([uniform(-r_spread, r_spread) for _ in range(3)]) for _ in range(self.n)]
        self.v = [np.array([uniform(-v_spread, v_spread) for _ in range(3)]) for _ in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        for i in range(self.n):
            self.q[i] /= np.linalg.norm(self.q[i])
        self.line = [[] for _ in range(self.n)]

        # Прорисовка ножек
        self.legs_x = 0.85
        self.legs_z = 0.7

class Objects:
    def __init__(self, n_c: int = 1, n_f: int = 5, model_c: str = '1U', if_any_print: bool = True):
        """Класс содержит информацию о n_c кубсатах и n_f фемтосатах. \n
        Размер n_c кубсатов определяется моделью model_c"""
        self.dt = 0.1
        self.n_a = 1
        self.n_f = n_f
        self.model_c = model_c

        # Классы
        self.p = PhysicModel()

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
