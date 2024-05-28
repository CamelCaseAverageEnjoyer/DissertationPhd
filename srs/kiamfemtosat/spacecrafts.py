"""Функции, связанные с архитектурой КА"""
from srs.kiamfemtosat.my_math import *
from srs.kiamfemtosat.config import *
from srs.kiamfemtosat.cosmetic import my_print

# >>>>>>>>>>>> Небольшие функции <<<<<<<<<<<<
def local_dipole(r, ind: str = 'x') -> float:
    """Возвращает диаграмму направленности одной полуволновой антенны (бублик). Костыль: возвращается >= 0
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
    tmp = np.cos(cos_a * np.pi / 2) / sin_a + DISTORTION * cos_a ** 2 + DISTORTION * aside
    if tmp < 0 and IF_ANY_PRINT:
        my_print(f"Внимание! Отрицательное усиление! G={tmp} (cos={cos_a}, sin={sin_a})", color="r")
    return clip(tmp, 0, 1e10)

def get_gain(obj: any, r: Union[float, np.ndarray], mode3: bool = False) -> list:
    """Внимание! Всё переделано - теперь возвращается только список для повышения градуса полиморфизма,
    интуитивизма, индуизма, культуризма, конституционализма, шовинизма, каннибализма"""
    # Памятка: GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
    r1 = r / np.linalg.norm(r)
    if obj.gain_mode == GAIN_MODES[1]:
        return [np.linalg.norm([r1[0] * 1, r1[1] * 0.7, r1[2] * 0.8])]
    if obj.gain_mode == GAIN_MODES[2]:
        return [local_dipole(r1, 'z')]
    if obj.gain_mode == GAIN_MODES[3] or (mode3 and obj.gain_mode == GAIN_MODES[4]):
        return [local_dipole(r1, 'x') + local_dipole(r1, 'y')]
    if obj.gain_mode == GAIN_MODES[4] and not mode3:
        return [local_dipole(r1, 'x'), local_dipole(r1, 'y')]
    if obj.gain_mode == GAIN_MODES[5] and not mode3:
        return [local_dipole(r1, 'x'), local_dipole(r1, 'y'), local_dipole(r1, 'z')]
    return [1]

# >>>>>>>>>>>> Классы аппаратов <<<<<<<<<<<<
class FemtoSat:
    def __init__(self, w_orb: float, n: int = 10, r_spread: float = 1e2, v_spread: float = 1e-1,
                 start_navigation: str = 'near', start_navigation_tolerance: float = 0.9, w_spread: float = 1e-5):
        """Класс содержит информацию об n фемтосатах.\n
        Все величны представлены в СИ."""
        # Общие параметры
        self.n = n
        self.mass = 0.01
        self.size = [0.03, 0.03]
        self.power_signal_full = 0.01
        self.length_signal_full = 0.001
        self.gain_mode = GAIN_MODES[0]

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-r_spread, r_spread, 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v_spread, v_spread, 3) for _ in range(self.n)]
        self.w_orf = [np.random.uniform(-w_spread, w_spread, 3) for _ in range(self.n)]
        self.q, self.q_ = [[np.random.uniform(-1, 1, 4) for _ in range(self.n)] for _ in range(2)]
        self.c_hkw, self.line, self.line_kalman, self.line_difference, self.attitude_difference, self.spin_difference, \
            self.z_difference = [[[] for _ in range(self.n)] for _ in range(7)]
        for i in range(self.n):
            if SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
            self.q[i] /= np.linalg.norm(self.q[i])
            self.q_[i] /= np.linalg.norm(self.q_[i])

        # Индивидуальные параметры режимов работы
        self.operating_mode = [OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [CHIPSAT_OPERATING_MODE for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_ = \
            [[[[] for _ in range(self.n)] for _ in range(self.n)] for _ in range(4)]
        prm_poor = [np.append(np.append(np.append(np.random.uniform(-r_spread, r_spread, 3), self.q_[i]),
                                        np.random.uniform(-v_spread, v_spread, 3)),
                              np.random.uniform(-w_spread, w_spread, 3)) for i in range(self.n)]
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i]), self.v_orf[i]), self.w_orf[i])
                    for i in range(self.n)]
        start_navigation_tolerance = 1 if start_navigation == NAVIGATIONS[0] else start_navigation_tolerance
        start_navigation_tolerance = 0 if start_navigation == NAVIGATIONS[2] else start_navigation_tolerance
        self.rv_orf_calc = [prm_good[i] * start_navigation_tolerance +
                            prm_poor[i] * (1 - start_navigation_tolerance) for i in range(self.n)]

class CubeSat:
    def __init__(self, w_orb: float, n_f: int, n: int = 1, model: str = '1U', w_spread: float = 1e-4,
                 r_spread: float = R_V_CubeSat_SPREAD[0], v_spread: float = R_V_CubeSat_SPREAD[1]):
        """Класс содержит информацию об n кубсатах модели model_c = 1U/1.5U/2U/3U/6U/12U.\n
        Все величны представлены в СИ."""
        # Предопределённые параметры
        models = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        masses = [2., 3., 4., 6., 12., 24.]
        mass_center_errors = [[0.02, 0.02, 0.02], [0.02, 0.02, 0.03], [0.02, 0.02, 0.045],
                              [0.02, 0.02, 0.07], [4.5, 2., 7.], [4.5, 4.5, 7.]]
        sizes = [[0.1, 0.1, 0.1135], [0.1, 0.1, 0.1702], [0.1, 0.1, 0.227],
                 [0.1, 0.1, 0.3405], [0.2263, 0.1, 0.366], [0.2263, 0.2263, 0.366]]
        if model not in models:
            raise ValueError(f"Модель кубсата [{model}] должна быть среди {models}")

        # Общие параметры
        self.n = n
        self.model = model
        self.model_number = models.index(model)
        self.mass = masses[self.model_number]
        self.mass_center_error = mass_center_errors[self.model_number]
        self.r_mass_center = np.array([np.random.uniform(-i, i) for i in self.mass_center_error])
        self.size = sizes[self.model_number]
        self.gain_mode = GAIN_MODES[0]

        # Индивидуальные параметры движения
        self.r_orf = [np.random.uniform(-r_spread, r_spread, 3) for _ in range(self.n)]
        self.v_orf = [np.random.uniform(-v_spread, v_spread, 3) for _ in range(self.n)]
        self.w_orf = [np.random.uniform(-w_spread, w_spread, 3) for _ in range(self.n)]
        self.q = [np.array([np.random.uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        self.c_hkw, self.line = [[[] for _ in range(self.n)] for _ in range(2)]
        for i in range(self.n):
            if SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
            self.q[i] /= np.linalg.norm(self.q[i])

        # Индивидуальные параметры режимов работы
        self.operating_mode = [OPERATING_MODES[0] for _ in range(self.n)]
        self.operating_modes = [OPERATING_MODES_CHANGE[0] for _ in range(self.n)]

        # Индивидуальные параметры управления
        self.m_self, self.b_env = [[np.zeros(3) for _ in range(self.n)] for _ in range(2)]

        # Индивидуальные параметры измерений
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_, self.kalm_dist = \
            [[[[] for _ in range(n_f)] for _ in range(self.n)] for _ in range(5)]

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

