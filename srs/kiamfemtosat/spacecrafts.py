"""Функции, связанные с архитектурой КА"""
from srs.kiamfemtosat.my_math import *
from srs.kiamfemtosat.config import Variables
from srs.kiamfemtosat.cosmetic import my_print

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

def get_gain(v: Variables, obj: any, r: Union[float, np.ndarray], mode3: bool = False) -> list:
    """Внимание! Всё переделано - теперь возвращается только список для повышения градуса полиморфизма,
    интуитивизма, индуизма, культуризма, конституционализма, шовинизма, каннибализма
    :param v: Объект класса Variables
    :param obj: Переменная класса FemtoSat или CubeSat (any потому что не хочу ниже писать)
    :param r: Направление сигнала в СК антенны
    :param mode3: Специальный флаг для возврата списка длины 1"""
    # Памятка: GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
    e = r / np.linalg.norm(r)
    if obj.gain_mode == v.GAIN_MODES[1]:
        return [np.linalg.norm([e[0] * 1, e[1] * 0.7, e[2] * 0.8])]
    if obj.gain_mode == v.GAIN_MODES[2]:
        return [local_dipole(v, e, 'z')]
    if obj.gain_mode == v.GAIN_MODES[3] or ((mode3 or v.MULTI_ANTENNA_USE) and obj.gain_mode == v.GAIN_MODES[4]):
        return [local_dipole(v, e, 'x') + local_dipole(v, e, 'y')]
    if obj.gain_mode == v.GAIN_MODES[4] and not (mode3 or v.MULTI_ANTENNA_USE):
        return [local_dipole(v, e, 'x'), local_dipole(v, e, 'y')]
    if (mode3 or v.MULTI_ANTENNA_USE) and obj.gain_mode == v.GAIN_MODES[5]:
        return [local_dipole(v, e, 'x') + local_dipole(v, e, 'y') + local_dipole(v, e, 'z')]
    if obj.gain_mode == v.GAIN_MODES[5] and not (mode3 or v.MULTI_ANTENNA_USE):
        return [local_dipole(v, e, 'x'), local_dipole(v, e, 'y'), local_dipole(v, e, 'z')]
    return [1]

# >>>>>>>>>>>> Классы аппаратов <<<<<<<<<<<<
class Anchor:
    def __init__(self, v: Variables):
        """Класс фантомного КА, движущегося по орбите с нулевым разбросом скоростей и положений"""
        from srs.kiamfemtosat.dynamics import get_matrices, o_i, get_c_hkw

        # Общие параметры
        self.name = "Anchor"
        self.n = 1
        self.mass = 1e10
        self.size = [1., 1., 1.]

        # Индивидуальные параметры движения
        self.w_irf = [np.zeros(3) for _ in range(self.n)]
        self.q = [np.array([1, 0, 0, 0]) for _ in range(self.n)]
        # self.q = [[np.array([1, 0, 0, 0]) for _ in range(self.n)] for _ in range(1)]
        # print(self.q)
        self.r_orf = [np.zeros(3) for _ in range(self.n)]
        self.v_orf = [np.zeros(3) for _ in range(self.n)]
        U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=0, first_init=True)
        self.r_irf = [o_i(v=v, a=self.r_orf[0], U=U, vec_type='r')]
        self.v_irf = [o_i(v=v, a=self.v_orf[0], U=U, vec_type='v')]
        self.line_orf, self.line_irf, self.line_kalman, self.line_difference, self.attitude_difference, \
            self.spin_difference, self.z_difference = [[[] for _ in range(self.n)] for _ in range(7)]
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры измерений
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i]), self.v_orf[i]), self.w_irf[i])
                    for i in range(self.n)]
        self.rv_orf_calc = [prm_good[i] for i in range(self.n)]

class FemtoSat:
    def __init__(self, v: Variables):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        from srs.kiamfemtosat.dynamics import get_matrices, o_i, get_c_hkw
        if v.CHIPSAT_AMOUNT < 0:
            raise ValueError(f"Количество чипсатов {v.CHIPSAT_AMOUNT} должно быть не меньше 0!")

        # Общие параметры
        self.name = "FemtoSat"
        self.n = v.CHIPSAT_AMOUNT
        self.mass = 0.01
        self.size = [0.03, 0.03]
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001
        self.gain_mode = v.GAIN_MODEL_F

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[0], v.RVW_ChipSat_SPREAD[0], 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[1], v.RVW_ChipSat_SPREAD[1], 3) for _ in range(self.n)]
        self.r_irf = [np.zeros(3) for _ in range(self.n)]  # Инициализируется ниже
        self.v_irf = [np.zeros(3) for _ in range(self.n)]
        self.w_irf = [np.random.uniform(-v.RVW_ChipSat_SPREAD[2], v.RVW_ChipSat_SPREAD[2], 3) for _ in range(self.n)]
        self.q, self.q_ = [[np.random.uniform(-1, 1, 4) for _ in range(self.n)] for _ in range(2)]
        self.line_orf, self.line_irf, self.line_kalman, self.line_difference, self.attitude_difference, \
            self.spin_difference, self.z_difference = [[[] for _ in range(self.n)] for _ in range(7)]
        for i in range(self.n):
            if v.SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * v.W_ORB
            self.q[i] /= np.linalg.norm(self.q[i])
            self.q_[i] /= np.linalg.norm(self.q_[i])
            U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=i)
            self.r_irf[i] = o_i(v=v, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(v=v, a=self.v_orf[i], U=U, vec_type='v')
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [v.OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [v.CHIPSAT_OPERATING_MODE for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_ = \
            [[[[] for _ in range(self.n)] for _ in range(self.n)] for _ in range(4)]
        prm_poor = [np.append(
            np.append(np.append(np.random.uniform(-v.RVW_ChipSat_SPREAD[0], v.RVW_ChipSat_SPREAD[0], 3), self.q_[i]),
                      np.random.uniform(-v.RVW_ChipSat_SPREAD[1], v.RVW_ChipSat_SPREAD[1], 3)),
            np.random.uniform(-v.RVW_ChipSat_SPREAD[2], v.RVW_ChipSat_SPREAD[2], 3)) for i in range(self.n)]
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i]), self.v_orf[i]), self.w_irf[i])
                    for i in range(self.n)]
        start_navigation_tolerance = 1 if v.START_NAVIGATION == v.NAVIGATIONS[0] else v.START_NAVIGATION_TOLERANCE
        start_navigation_tolerance = 0 if v.START_NAVIGATION == v.NAVIGATIONS[2] else start_navigation_tolerance
        self.rv_orf_calc = [prm_good[i] * start_navigation_tolerance +
                            prm_poor[i] * (1 - start_navigation_tolerance) for i in range(self.n)]

class CubeSat:
    def __init__(self, v: Variables):
        """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.\n
        Все величны представлены в СИ."""
        from srs.kiamfemtosat.dynamics import get_matrices, o_i, get_c_hkw
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
        self.v_irf = [np.zeros(3) for _ in range(self.n)]
        self.w_irf = [np.random.uniform(-v.RVW_CubeSat_SPREAD[2], v.RVW_CubeSat_SPREAD[2], 3) for _ in range(self.n)]
        self.q = [np.array([np.random.uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        self.line_orf, self.line_irf = [[[] for _ in range(self.n)] for _ in range(2)]
        for i in range(self.n):
            if v.SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * v.W_ORB
            self.q[i] /= np.linalg.norm(self.q[i])
            U, _, _, _ = get_matrices(v=v, t=0, obj=self, n=i)
            self.r_irf[i] = o_i(v=v, a=self.r_orf[i], U=U, vec_type='r')
            self.v_irf[i] = o_i(v=v, a=self.v_orf[i], U=U, vec_type='v')
        self.c_hkw = [get_c_hkw(self.r_orf[i], self.v_orf[i], v.W_ORB) for i in range(self.n)]

        # Индивидуальные параметры режимов работы
        self.operating_mode = [v.OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [v.OPERATING_MODES_CHANGE[0] for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_, self.kalm_dist = \
            [[[[] for _ in range(v.CHIPSAT_AMOUNT)] for _ in range(self.n)] for _ in range(5)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007
