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
    def __init__(self, v: Variables):
        """Класс фантомного КА, движущегося по орбите с нулевым разбросом скоростей и положений"""
        from dynamics import get_matrices, o_i, get_c_hkw

        # Общие параметры
        self.name = "No exist"
        self.n = 1
        self.mass = 1e10
        self.size = [1., 1., 1.]

        # Индивидуальные параметры движения
        self.w_irf = [np.zeros(3) for _ in range(self.n)]  # Убрать потом !!!!!!!!!!!!!!
        self.w_orf = [np.zeros(3) for _ in range(self.n)]  # Убрать потом !!!!!!!!!!!!!!
        self.q = [np.array([1, 0, 0, 0]) for _ in range(self.n)]
        self.r_orf = [np.zeros(3) for _ in range(self.n)]
        self.v_orf = [np.zeros(3) for _ in range(self.n)]
        U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=0, first_init=True)
        self.r_irf = [o_i(v=v, a=self.r_orf[0], U=U, vec_type='r')]
        self.v_irf = [o_i(v=v, a=self.v_orf[0], U=U, vec_type='v')]
        self.line_orf, self.line_irf, self.line_kalman, self.line_difference, self.attitude_difference, \
            self.spin_difference = [[[] for _ in range(self.n)] for _ in range(6)]
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры измерений
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i][1:4]), self.v_orf[i]), self.w_irf[i])
                    for i in range(self.n)]
        self.rv_orf_calc = [prm_good[i] for i in range(self.n)]


class Anchor(Apparatus):
    def __init__(self, v: Variables):
        """Класс фантомного КА, движущегося по орбите с нулевым разбросом скоростей и положений"""
        super().__init__(v)
        self.name = "Anchor"

class FemtoSat(Apparatus):
    def __init__(self, v: Variables):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        from dynamics import get_matrices, o_i, get_c_hkw
        if v.CHIPSAT_AMOUNT < 0:
            raise ValueError(f"Количество чипсатов {v.CHIPSAT_AMOUNT} должно быть не меньше 0!")

        # Предопределённые параметры
        masses = [0.01, 0.1]
        mass_center_errors = [[0.001, -0.001], [0.005, 0.003]]
        sizes = [[0.03, 0.03], [0.4, 0.15]]

        # Общие параметры
        self.name = "FemtoSat"
        self.n = v.CHIPSAT_AMOUNT
        self.model = v.CHIPSAT_MODEL
        self.model_number = v.CHIPSAT_MODELS.index(self.model)
        self.mass = masses[self.model_number]
        self.mass_center_error = mass_center_errors[self.model_number]
        self.size = sizes[self.model_number]
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001
        self.gain_mode = v.GAIN_MODEL_F

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[0], v.RVW_ChipSat_SPREAD[0], 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[1], v.RVW_ChipSat_SPREAD[1], 3) for _ in range(self.n)]
        self.r_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.v_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.w_orf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[2], v.RVW_ChipSat_SPREAD[2], 3) for _ in range(self.n)]
        self.w_orf_ = [np.random.uniform(-v.RVW_ChipSat_SPREAD[2], v.RVW_ChipSat_SPREAD[2], 3) for _ in range(self.n)]
        self.w_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.w_irf_ = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.q, self.q_ = [[np.random.uniform(-1, 1, 4) for _ in range(self.n)] for _ in range(2)]
        self.line_orf, self.line_irf, self.line_kalman, self.line_difference, self.attitude_difference, \
            self.spin_difference = [[[] for _ in range(self.n)] for _ in range(6)]
        for i in range(self.n):
            if v.SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * v.W_ORB
            self.q[i] /= np.linalg.norm(self.q[i])
            self.q_[i] /= np.linalg.norm(self.q_[i])
            U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=i)
            self.r_irf[i] = o_i(v=v, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(v=v, a=self.v_orf[i], U=U, vec_type='v')
            self.w_irf[i] = o_i(a=self.w_orf[i], v=v, U=U, vec_type='w')
            self.w_irf_[i] = o_i(a=self.w_orf_[i], v=v, U=U, vec_type='w')
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [v.OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [v.CHIPSAT_OPERATING_MODE for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.dist_estimate, self.dist_based_measures = \
            [[[[] for _ in range(self.n)] for _ in range(self.n)] for _ in range(4)]
        prm_poor = [np.append(
            np.append(np.append(np.random.uniform(-v.RVW_ChipSat_SPREAD[0], v.RVW_ChipSat_SPREAD[0], 3),
                                self.q_[i][1:4]), self.w_irf_[i]),
            np.random.uniform(-v.RVW_ChipSat_SPREAD[2], v.RVW_ChipSat_SPREAD[2], 3)) for i in range(self.n)]
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i][1:4]), self.v_orf[i]), self.w_orf[i])
                    for i in range(self.n)]
        start_navigation_tolerance = 1 if v.START_NAVIGATION == v.NAVIGATIONS[0] else v.START_NAVIGATION_TOLERANCE
        start_navigation_tolerance = 0 if v.START_NAVIGATION == v.NAVIGATIONS[2] else start_navigation_tolerance
        self.rv_orf_calc = [prm_good[i] * start_navigation_tolerance +
                            prm_poor[i] * (1 - start_navigation_tolerance) for i in range(self.n)]

class CubeSat(Apparatus):
    """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.
    Все величны представлены в СИ."""
    def __init__(self, v: Variables):
        from dynamics import get_matrices, o_i, get_c_hkw
        if v.CUBESAT_AMOUNT < 1:
            raise ValueError(f"Количество чипсатов {v.CUBESAT_AMOUNT} должно быть не меньше 1!")

        # Предопределённые параметры
        masses = [2., 3., 4., 6., 12., 24.]
        mass_center_errors = [[0.02, 0.02, 0.02], [0.02, 0.02, 0.03], [0.02, 0.02, 0.045],
                              [0.02, 0.02, 0.07], [4.5, 2., 7.], [4.5, 4.5, 7.]]
        sizes = [[0.1, 0.1, 0.1135], [0.1, 0.1, 0.1702], [0.1, 0.1, 0.227],
                 [0.1, 0.1, 0.3405], [0.2263, 0.1, 0.366], [0.2263, 0.2263, 0.366]]

        # Общие параметры
        self.name = "CubeSat"
        self.n = v.CUBESAT_AMOUNT
        self.model = v.CUBESAT_MODEL
        self.model_number = v.CUBESAT_MODELS.index(self.model)
        self.mass = masses[self.model_number]
        self.mass_center_error = mass_center_errors[self.model_number]
        self.r_mass_center = np.array([np.random.uniform(-i, i) for i in self.mass_center_error])
        self.size = sizes[self.model_number]
        self.gain_mode = v.GAIN_MODEL_C

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[0], v.RVW_CubeSat_SPREAD[0], 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[1], v.RVW_CubeSat_SPREAD[1], 3) for _ in range(self.n)]
        self.r_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.v_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.w_orf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[2], v.RVW_CubeSat_SPREAD[2], 3) for _ in range(self.n)]
        self.w_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.q = [np.array([np.random.uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        self.line_orf, self.line_irf = [[[] for _ in range(self.n)] for _ in range(2)]
        for i in range(self.n):
            if v.SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * v.W_ORB
            self.q[i] /= np.linalg.norm(self.q[i])
            U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=i)
            self.r_irf[i] = o_i(v=v, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(v=v, a=self.v_orf[i], U=U, vec_type='v')
            self.w_irf[i] = o_i(a=self.w_orf[i], v=v, U=U, vec_type='w')
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [v.OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [v.OPERATING_MODES_CHANGE[0] for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.dist_estimate, self.dist_based_measures, self.kalm_dist = \
            [[[[] for _ in range(v.CHIPSAT_AMOUNT)] for _ in range(self.n)] for _ in range(5)]
        self.z_difference = [[[] for _ in range(v.CHIPSAT_AMOUNT)] for _ in range(self.n)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007
