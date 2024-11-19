"""Функции, связанные с архитектурой КА"""
from my_math import *
from config import Variables
from cosmetic import my_print

# >>>>>>>>>>>> Диаграмма направленности антенн связи <<<<<<<<<<<<
def local_dipole(v: Variables, r: Union[list, np.ndarray], ind: str = 'x') -> float:
    """Возвращает диаграмму направленности одной полуволновой антенны (бублик). Костыль: возвращается >= 0
    :param v: Объект класса Variables
    :param r: Радиус-вектор направления антенны
    :param ind: Координата направления антенны"""
    if ind not in "xyz":
        raise ValueError(f"Координата «{ind}» должна быть среди: [x, y, z]")
    r_m = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    cos_a = (r[0] * int(ind == 'x') + r[1] * int(ind == 'y') + r[2] * int(ind == 'z')) / r_m
    sin_a = (np.sqrt(r[1] ** 2 + r[2] ** 2) * int(ind == 'x') +
             np.sqrt(r[0] ** 2 + r[2] ** 2) * int(ind == 'y') +
             np.sqrt(r[0] ** 2 + r[1] ** 2) * int(ind == 'z')) / r_m
    aside = ((r[1] + r[2]) * int(ind == 'x') +
             (r[0] + r[2]) * int(ind == 'y') +
             r[2] * int(ind == 'z')) / r_m
    tmp = np.cos(cos_a * np.pi / 2) / sin_a + v.DISTORTION * cos_a ** 2 + v.DISTORTION * aside
    my_print(f"Внимание! Отрицательное усиление! G={tmp} (cos={cos_a}, sin={sin_a})", color="r",
             if_print=tmp < 0 and v.IF_ANY_PRINT)
    return clip(tmp, 0, 1e10)

def get_gain(v: Variables, obj: any, r: Union[float, np.ndarray], if_take: bool = False, if_send: bool = False) -> list:
    """Внимание! Всё переделано - теперь возвращается только список для повышения градуса полиморфизма,
    интуитивизма, индуизма, культуризма, конституционализма, шовинизма, каннибализма
    :param v: Объект класса Variables
    :param obj: Переменная класса FemtoSat или CubeSat (any потому что не хочу ниже писать)
    :param r: Направление сигнала в СК антенны
    :param if_take: Флаг на принятый сигнал
    :param if_send: Флаг на посланный сигнал"""
    # Памятка: GAIN_MODES = ['isotropic', '1 antenna', '2 antennas', '3 antennas', 'ellipsoid']
    e = r / np.linalg.norm(r)
    if obj.gain_mode == v.GAIN_MODES[1]:
        return [local_dipole(v, e, 'x')]
    if obj.gain_mode == v.GAIN_MODES[2]:
        if (if_take and v.MULTI_ANTENNA_TAKE) or (if_send and v.MULTI_ANTENNA_SEND):
            return [local_dipole(v, e, 'x'), local_dipole(v, e, 'y')]
        return [local_dipole(v, e, 'x') + local_dipole(v, e, 'y')]
    if obj.gain_mode == v.GAIN_MODES[3]:
        if (if_take and v.MULTI_ANTENNA_TAKE) or (if_send and v.MULTI_ANTENNA_SEND):
            return [local_dipole(v, e, 'x'), local_dipole(v, e, 'y'), local_dipole(v, e, 'z')]
        return [local_dipole(v, e, 'x') + local_dipole(v, e, 'y') + local_dipole(v, e, 'z')]
    if obj.gain_mode == v.GAIN_MODES[4]:
        return [np.linalg.norm([e[0] * 1, e[1] * 0.7, e[2] * 0.8])]
    return [1]


# >>>>>>>>>>>> Классы аппаратов <<<<<<<<<<<<
class Apparatus:
    def __init__(self, v: Variables, n: int):
        """Пустой класс КА"""
        from dynamics import get_matrices, o_i, get_c_hkw

        # Общие параметры
        self.name = "No exist"
        self.n = n
        self.mass = 1e10
        self.size = [1., 1., 1.]
        self.J = None

        # Индивидуальные параметры движения
        self.w_irf = [np.zeros(3) for _ in range(self.n)]
        self.w_orf = [np.zeros(3) for _ in range(self.n)]
        self.q = [np.quaternion(1, 0, 0, 0) for _ in range(self.n)]
        self.r_orf = [np.zeros(3) for _ in range(self.n)]
        self.v_orf = [np.zeros(3) for _ in range(self.n)]
        U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=0, first_init=True)
        self.r_irf = [o_i(v=v, a=self.r_orf[0], U=U, vec_type='r')]
        self.v_irf = [o_i(v=v, a=self.v_orf[0], U=U, vec_type='v')]
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [v.OPERATING_MODES[0] for _ in range(self.n)]

        # Индивидуальные параметры измерений
        '''
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i][1:4]), self.v_orf[i]), self.w_irf[i])
                    for i in range(self.n)]
        self.rv_orf_calc = [prm_good[i] for i in range(self.n)]'''

    def update_c(self, v):
        from dynamics import get_c_hkw
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

    def update_irf_rv(self, v: Variables, t: float = 0):
        from dynamics import o_i, get_matrices
        for i in range(self.n):
            U, _, _, _ = get_matrices(v=v, t=t, obj=self, n=i)
            self.r_irf[i] = o_i(v=v, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(v=v, a=self.v_orf[i], U=U, vec_type='v')

    def update_irf_w(self, v: Variables, t: float = 0, w_irf: list = None, w_orf: list = None):
        from dynamics import o_i, get_matrices
        w_irf = self.w_irf if w_irf is None else w_irf
        w_orf = self.w_orf if w_orf is None else w_orf
        for i in range(self.n):
            U, _, _, _ = get_matrices(v=v, t=t, obj=self, n=i)
            w_irf[i] = o_i(a=w_orf[i], v=v, U=U, vec_type='w')

    def init_correct_q_v(self, v: Variables, q: list = None):
        q = self.q if q is None else q
        for i in range(self.n):
            if v.SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * v.W_ORB
            q[i] = q[i].normalized()


class Anchor(Apparatus):
    def __init__(self, v: Variables):
        """Класс фантомного КА, движущегося по орбите с нулевым разбросом скоростей и положений"""
        super().__init__(v=v, n=1)
        self.name = "Anchor"

class CubeSat(Apparatus):
    """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.
    Все величны представлены в СИ."""
    def __init__(self, v: Variables):
        super().__init__(v=v, n=v.CUBESAT_AMOUNT)

        # Предопределённые параметры
        cubesat_property = {'1U': {'mass': 2.,
                                   'mass_center_error': [0.02, 0.02, 0.02],
                                   'dims': [0.1, 0.1, 0.1135]},
                            '1.5U': {'mass': 3.,
                                     'mass_center_error': [0.02, 0.02, 0.03],
                                     'dims': [0.1, 0.1, 0.1702]},
                            '2U': {'mass': 4.,
                                   'mass_center_error': [0.02, 0.02, 0.045],
                                   'dims': [0.1, 0.1, 0.227]},
                            '3U': {'mass': 6.,
                                   'mass_center_error': [0.02, 0.02, 0.07],
                                   'dims': [0.1, 0.1, 0.3405]},
                            '6U': {'mass': 12.,
                                   'mass_center_error': [4.5, 2., 7.],
                                   'dims': [0.2263, 0.1, 0.366]},
                            '12U': {'mass': 24.,
                                    'mass_center_error': [4.5, 4.5, 7.],
                                    'dims': [0.2263, 0.2263, 0.366]}}

        # Общие параметры
        self.name = "CubeSat"
        self.gain_mode = v.GAIN_MODEL_C
        self.mass = cubesat_property[v.CUBESAT_MODEL]['mass']
        self.size = cubesat_property[v.CUBESAT_MODEL]['dims']
        self.mass_center_error = cubesat_property[v.CUBESAT_MODEL]['mass_center_error']
        self.r_mass_center = np.array([np.random.uniform(-i, i) for i in self.mass_center_error])
        # Пока что J диагонален
        self.J = np.diag([self.size[1]**2 + self.size[2]**2,
                          self.size[0]**2 + self.size[2]**2,
                          self.size[0]**2 + self.size[1]**2]) * self.mass / 12

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[0], v.RVW_CubeSat_SPREAD[0], 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[1], v.RVW_CubeSat_SPREAD[1], 3) for _ in range(self.n)]
        self.w_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[2], v.RVW_CubeSat_SPREAD[2], 3) for _ in range(self.n)]
        # Инициализируется автоматически
        # self.q = [np.quaternion(1, 0, 0, 0) for _ in range(self.n)]
        self.q = [np.quaternion(*np.random.uniform(-1, 1, 4)) for _ in range(self.n)]
        self.init_correct_q_v(v=v)
        self.r_irf, self.v_irf, self.w_irf = [[np.zeros(3) for _ in range(self.n)] for _ in range(3)]
        self.update_irf_rv(v=v, t=0)
        self.update_irf_w(v=v, t=0)
        self.update_c(v=v)

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

class FemtoSat(Apparatus):
    def __init__(self, v: Variables, c: CubeSat):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        super().__init__(v=v, n=v.CHIPSAT_AMOUNT)

        # Предопределённые параметры
        chipsat_property = {'KickSat': {'mass': 0.01,
                                        'mass_center_error': [0.001, -0.001],
                                        'dims': [0.03, 0.03]},
                            '1.Трисат': {'mass': 0.1,
                                         'mass_center_error': [0.005, 0.003],
                                         'dims': [0.4, 0.15]}}

        # Общие параметры
        self.name = "FemtoSat"
        self.gain_mode = v.GAIN_MODEL_F
        self.mass = chipsat_property[v.CHIPSAT_MODEL]['mass']
        self.mass_center_error = chipsat_property[v.CHIPSAT_MODEL]['mass_center_error']
        self.size = chipsat_property[v.CHIPSAT_MODEL]['dims']
        # Пока что J диагонален
        self.J = np.diag([self.size[1]**2, self.size[0]**2, self.size[0]**2 + self.size[1]**2]) * self.mass / 12
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001

        def spread(_i: int):
            return np.random.uniform(-v.RVW_ChipSat_SPREAD[_i], v.RVW_ChipSat_SPREAD[_i], 3)

        # Индивидуальные параметры движения
        self.deploy(v=v, c=c, i_c=0, spread=spread)
        self.w_orf_ = [spread(2) for _ in range(self.n)]
        # Инициализируется автоматически
        self.r_irf, self.v_irf, self.w_irf, self.w_irf_ = [[np.zeros(3) for _ in range(self.n)] for _ in range(4)]
        self.q, self.q_ = [[np.quaternion(*np.random.uniform(-1, 1, 4)) for _ in range(self.n)] for _ in range(2)]
        self.init_correct_q_v(v=v)
        self.init_correct_q_v(v=v, q=self.q_)
        self.update_irf_rv(v=v, t=0)
        self.update_irf_w(v=v, t=0)
        self.update_irf_w(v=v, t=0, w_irf=self.w_irf_, w_orf=self.w_orf_)
        self.update_c(v=v)

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        tol = 1 if v.START_NAVIGATION == v.NAVIGATIONS[0] else v.START_NAVIGATION_TOLERANCE
        tol = 0 if v.START_NAVIGATION == v.NAVIGATIONS[2] else tol

        # Новый формат
        self.apriori_params = {'r orf': [self.r_orf[i] * tol + spread(0) * (1 - tol) for i in range(self.n)],
                               'v orf': [self.v_orf[i] * tol + spread(1) * (1 - tol) for i in range(self.n)],
                               'w irf': [self.w_irf[i] * tol + self.w_irf_[i] * (1 - tol) for i in range(self.n)],
                               'q-3 irf': [self.q[i].vec * tol + self.q_[i].vec * (1 - tol) for i in range(self.n)]}

    def deploy(self, v: Variables, c: CubeSat, i_c: int, spread) -> None:
        """
        Функция отделения задаёт начальные условия для дочерних КА из материнских КА
        :param v: объект Variables
        :param c: объект CubeSat
        :param i_c: id-номер материнского КА, от которого отделяются дочерние КА
        :param spread: функция np.random.uniform
        :return: {'r orf': ..., 'v orf': ..., 'q-3 irf': ..., 'w irf': ...}, где значения - list of np.ndarray
        """
        if v.DEPLOYMENT == v.DEPLOYMENTS[0]:  # No
            self.r_orf = [spread(0) for _ in range(self.n)]
            self.v_orf = [spread(1) for _ in range(self.n)]
            self.w_orf = [spread(2) for _ in range(self.n)]
        elif v.DEPLOYMENT == v.DEPLOYMENTS[1]:  # Specific
            r_before = c.r_orf[i_c]
            v_before = c.v_orf[i_c]
            dv = 1e-2
            v_deploy = v.RVW_ChipSat_SPREAD[1]
            self.r_orf = [r_before.copy() for _ in range(self.n)]
            self.v_orf = [v_before + np.array([0, 0, v_deploy]) + np.random.uniform(-dv, dv, 3) for _ in range(self.n)]
            self.w_orf = [spread(2) for _ in range(self.n)]
