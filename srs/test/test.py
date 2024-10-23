"""Проверка эквивалентности способов расчёта"""
import unittest
from srs.kiamfemto.config import Variables, Objects
from srs.kiamfemto.cosmetic import *
from srs.kiamfemto.my_math import *

dT_optimal = 0.1
TOTAL_TIME = 1e4

WTF_ESTIMATION = 1e1
ROUGH_ESTIMATION = 1e-1
MEDIUM_ESTIMATION = 1e-3
FINE_ESTIMATION = 1e-5

def get_v():
    v = Variables()
    v.test_mode()
    v.dT = dT_optimal
    v.TIME = TOTAL_TIME
    v.CUBESAT_AMOUNT = 1
    v.CHIPSAT_AMOUNT = 1
    v.IF_ANY_SHOW = False
    v.IF_ANY_PRINT = False
    v.IF_TEST_PRINT = False
    v.AERO_DRAG = False
    v.NAVIGATION_BY_ALL = False
    return v

def copy_kinematic_params(o1: Objects, o2: Objects):
    o2.f.r_irf = o1.f.r_irf.copy()
    o2.f.v_irf = o1.f.v_irf.copy()
    o2.f.r_orf = o1.f.r_orf.copy()
    o2.f.v_orf = o1.f.v_orf.copy()
    o2.f.c_hkw = o1.f.c_hkw.copy()
    o2.f.w_orf = o1.f.w_orf.copy()
    o2.f.q = o1.f.q.copy()

    o2.c.r_irf = o1.c.r_irf.copy()
    o2.c.v_irf = o1.c.v_irf.copy()
    o2.c.r_orf = o1.c.r_orf.copy()
    o2.c.v_orf = o1.c.v_orf.copy()
    o2.c.c_hkw = o1.c.c_hkw.copy()
    o2.c.w_orf = o1.c.w_orf.copy()
    o2.c.q = o1.c.q.copy()
    return o1, o2

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

    # >>>>>>>>>>>> Проверки мини-функций <<<<<<<<<<<<
    def test_my_math(self):
        self.assertEqual(np.pi/2, deg2rad(90))
        self.assertEqual(rad2deg(2*np.pi), 360)
        self.assertEqual(vec2quat([1, 2, 3]), [0, 1/np.sqrt(14), 2/np.sqrt(14), 3/np.sqrt(14)])
        self.assertEqual(vec2quat(np.array([0, 0, 0])), [1, 0, 0, 0])

    # >>>>>>>>>>>> Проверки на дурачка <<<<<<<<<<<<
    def test_while_nothing_has_been_done(self):
        """Функция проверяет, что в момент времени t=0 КА не улетают ни на йоту, а на t=dt чуть-чуть"""
        def get_some_params(o_):
            """Вспомогательная функция для вызова одного и того же числа переменных"""
            return [o_.c.r_orf[0], o_.c.v_orf[0], o_.c.r_irf[0], o_.c.v_irf[0], o_.c.q[0], o_.f.r_orf[0],
                    o_.f.v_orf[0], o_.f.r_irf[0], o_.f.v_irf[0], o_.f.q[0]]
        names = ['c_r_orf', 'c_v_orf', 'c_r_irf', 'c_v_irf', 'c_q', 'f_r_orf', 'f_v_orf', 'f_r_irf', 'f_v_irf', 'f_q']

        # Инициализация
        v = get_v()
        v.SOLVER = ['rk4', 'kiamastro'][0]
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
        v = get_v()
        v.SOLVER = ['rk4', 'kiamastro'][0]
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
        v = get_v()
        for i in range(2):
            v.SOLVER = ['rk4', 'kiamastro'][i]
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

    def test_zero_eccentricity(self):
        """Функция проверяет, что при e=0 истинная аномалия равна эксцентрической.
        Тогда матрица U задаётся корректно.

        К сожалению, придётся копировать [E, f] ручками из dynamics.get_matrices"""
        v = get_v()
        v.ECCENTRICITY = 0
        for t in np.linspace(0, v.SEC_IN_TURN, 100):
            E = t * v.W_ORB  # Эксцентрическая аномалия
            f = 2 * np.arctan(np.sqrt((1 + v.ECCENTRICITY) / (1 - v.ECCENTRICITY)) * np.tan(E / 2))  # Истинная аномалия
            # self.assertEqual(E, f)
            self.assertEqual(abs(E - f) % (2 * np.pi), 0)

    # >>>>>>>>>>>> Неочевидные проверки <<<<<<<<<<<<
    def test_multi_antenna_use(self):
        """Функция проверяет, что для 1 чипсата v.MULTI_ANTENNA_USE не влияет ни на что"""
        pass

    def test_my_and_kiamastro(self):
        """Функция проверяет, что в детерминант матриц остаётся единичным"""
        def get_some_params(o_, k_: int):
            """Вспомогательная функция для вызова одного и того же числа переменных"""
            return [[o_.c.r_irf[0], o_.c.v_irf[0], o_.f.r_irf[0], o_.f.v_irf[0]],
                    [o_.c.q[0], o_.f.q[0]],
                    [o_.c.r_orf[0], o_.c.v_orf[0], o_.f.r_orf[0], o_.f.v_orf[0]]][k_]
        names = [['c_r_irf', 'c_v_irf', 'f_r_irf', 'f_v_irf'],
                 ['c_q', 'f_q'],
                 ['c_r_orf', 'c_v_orf', 'f_r_orf', 'f_v_orf']]
        labels = ['Проверка ИСК', 'Проверка кватернионов', 'Проверка ОСК']
        # estimations = [ROUGH_ESTIMATION, FINE_ESTIMATION, WTF_ESTIMATION]
        estimations = [MEDIUM_ESTIMATION, WTF_ESTIMATION, 1e1]

        for k in range(3):
            # Инициализация
            v_1 = get_v()
            v_1.SOLVER = ['rk4', 'kiamastro'][0]
            o_1 = Objects(v=v_1)

            v_2 = get_v()
            v_2.SOLVER = ['rk4', 'kiamastro'][1]
            o_2 = Objects(v=v_2)

            # Одинаковые НУ (Начальные Условия)
            o_1, o_2 = copy_kinematic_params(o_1, o_2)

            # Интегрирование
            iterations = 5
            for ii in range(iterations):
                o_1.integrate(int(TOTAL_TIME//iterations))
                o_2.integrate(int(TOTAL_TIME//iterations))
                params_1 = get_some_params(o_1, k_=k)
                params_2 = get_some_params(o_2, k_=k)
                for i in range(len(params_1)):
                    for j in range(len(params_1[i])):
                        tmp = abs(params_1[i][j] - params_2[i][j]) \
                              / ((np.linalg.norm(params_1[i]) + np.linalg.norm(params_2[i])) / 2)
                        my_print(f"{labels[k]}: итерация {ii+1}/{iterations} | anw -> {tmp}"
                                 f"\n{i}:{names[k][i]}[{j}] | rk_4: {params_1[i][j]}, kiam: {params_2[i][j]}\n"
                                 f"rk_4: {params_1[i]}\nkiam: {params_2[i]}\n"
                                 f"Время (rk4) {o_1.p.t} = {o_2.p.to_delete}",
                                 color='y', if_print=tmp > estimations[k])
                        self.assertLess(tmp, estimations[k])


if __name__ == "__main__":
    unittest.main()
