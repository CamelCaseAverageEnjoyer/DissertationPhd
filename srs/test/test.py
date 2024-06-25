"""Проверка эквивалентности способов расчёта"""
import unittest
import numpy as np
from srs.kiamfemtosat.main import Objects
from srs.kiamfemtosat.config import Variables
from srs.kiamfemtosat.cosmetic import *

dT_optimal = 1.
ROUGH_ESTIMATION = 1e-1
FINE_ESTIMATION = 1e-5


class MyTests(unittest.TestCase):
    """Методы:
    self.assertEqual
    self.assertTrue
    self.assertFalse

    s = 'hello world'
    self.assertEqual(s.split(), ['hello', 'world'])
    """
    # >>>>>>>>>>>> Проверки на рукописные ошибки <<<<<<<<<<<<
    def test_necessary_errors_1(self):
        """Функция проверяет, что код ругается на кривые переменные"""
        v = Variables()
        v.CUBESAT_AMOUNT = 1
        v.CHIPSAT_AMOUNT = -1
        with self.assertRaises(ValueError):
            _ = Objects(v=v)
        v.CUBESAT_AMOUNT = 0
        v.CHIPSAT_AMOUNT = 1
        with self.assertRaises(ValueError):
            _ = Objects(v=v)

    # >>>>>>>>>>>> Проверки на дурачка <<<<<<<<<<<<
    def test_while_nothing_has_been_done(self):
        """Функция проверяет, что в момент времени t=0 КА не улетают ни на йоту, а на t=dt чуть-чуть"""
        def get_some_params(o_):
            """Вспомогательная функция для вызова одного и того же числа переменных"""
            return [o_.c.r_orf[0], o_.c.v_orf[0], o_.c.r_irf[0], o_.c.v_irf[0], o_.c.q[0], o_.f.r_orf[0],
                    o_.f.v_orf[0], o_.f.r_irf[0], o_.f.v_irf[0], o_.f.q[0]]
        names = ['c_r_orf', 'c_v_orf', 'c_r_irf', 'c_v_irf', 'c_q', 'f_r_orf', 'f_v_orf', 'f_r_irf', 'f_v_irf', 'f_q']

        # Инициализация
        v = Variables()
        v.test_mode()
        v.dT = dT_optimal
        v.DYNAMIC_MODEL = ['rk4', 'kiamastro'][0]
        v.NAVIGATION_BY_ALL = False
        v.CUBESAT_AMOUNT = 1
        v.CHIPSAT_AMOUNT = 1
        o = Objects(v=v)
        params = get_some_params(o).copy()

        # Интегрирование на 0 секунд
        o.integrate(t=0)
        for i in range(len(params)):
            for j in range(len(params[i])):
                self.assertEqual(get_some_params(o)[i][j], params[i][j])

        # Интегрирование на 1 итерацию
        o.integrate(t=v.dT)
        for i in range(len(params)):
            for j in range(len(params[i])):
                tmp = abs((get_some_params(o)[i][j] - params[i][j]) /
                          (np.linalg.norm(get_some_params(o)[i]) + np.linalg.norm(params[i])) * 2)
                my_print(f"{i}:{names[i]}[{j}] | t=0: {params[i][j]}, t=dt: {get_some_params(o)[i][j]}\n"
                         f"t=0: {params[i]}\nt=dt: {get_some_params(o)[i]}",
                         color='y', if_print=tmp > ROUGH_ESTIMATION)
                self.assertLess(tmp, ROUGH_ESTIMATION)  # Порядок точности

    def test_limited_motion(self):
        """Функция проверяет, что в КА не улетают в бездну"""
        def get_some_params(o_):
            """Вспомогательная функция для вызова одного и того же числа переменных"""
            return [o_.c.r_irf[0], o_.c.v_irf[0], o_.c.q[0], o_.f.r_irf[0], o_.f.v_irf[0], o_.f.q[0]]

        # Инициализация
        v = Variables()
        v.test_mode()
        v.dT = dT_optimal
        v.DYNAMIC_MODEL = ['rk4', 'kiamastro'][0]
        v.NAVIGATION_BY_ALL = False
        v.CUBESAT_AMOUNT = 1
        v.CHIPSAT_AMOUNT = 1
        o = Objects(v=v)

        # Интегрирование
        o.integrate(t=1e4)
        params = get_some_params(o).copy()
        for i in range(len(params)):
            for j in range(len(params[i])):
                self.assertLess(abs(params[i][j]), 1e8)

    def test_limited_matrices(self):
        """Функция проверяет, что в детерминант матриц остаётся единичным"""
        from srs.kiamfemtosat.dynamics import get_matrices

        # Инициализация
        v = Variables()
        v.test_mode()
        v.dT = dT_optimal
        v.NAVIGATION_BY_ALL = False
        v.CUBESAT_AMOUNT = 1
        v.CHIPSAT_AMOUNT = 1
        for i in range(2):
            v.DYNAMIC_MODEL = ['rk4', 'kiamastro'][i]
            o = Objects(v=v)

            # Интегрирование
            for _ in range(10):
                o.integrate(t=1e3)
                for obj in [o.c, o.f]:
                    U, S, A, R_orb = get_matrices(v=v, t=o.p.t, obj=obj, n=0)
                    tmp = abs(1 - np.linalg.det(U))
                    self.assertLess(tmp, FINE_ESTIMATION)
                    tmp = abs(1 - np.linalg.det(S))
                    self.assertLess(tmp, FINE_ESTIMATION)
                    tmp = abs(1 - np.linalg.det(A))
                    self.assertLess(tmp, FINE_ESTIMATION)

    # >>>>>>>>>>>> Неочевидные проверки <<<<<<<<<<<<
    def test_my_and_kiamastro(self):
        """Функция проверяет, что в детерминант матриц остаётся единичным"""
        from srs.kiamfemtosat.dynamics import get_matrices

        # Инициализация
        v = Variables()
        v.test_mode()
        v.dT = dT_optimal
        v.NAVIGATION_BY_ALL = False
        v.CUBESAT_AMOUNT = 1
        v.CHIPSAT_AMOUNT = 1
        v.DYNAMIC_MODEL = ['rk4', 'kiamastro'][0]
        o_1 = Objects(v=v)

        '''# Интегрирование
        for _ in range(10):
            o.integrate(t=1e3)
            for obj in [o.c, o.f]:
                U, S, A, R_orb = get_matrices(v=v, t=o.p.t, obj=obj, n=0)
                tmp = abs(1 - np.linalg.det(U))
                self.assertLess(tmp, FINE_ESTIMATION)
                tmp = abs(1 - np.linalg.det(S))
                self.assertLess(tmp, FINE_ESTIMATION)
                tmp = abs(1 - np.linalg.det(A))
                self.assertLess(tmp, FINE_ESTIMATION)'''


if __name__ == "__main__":
    unittest.main()
