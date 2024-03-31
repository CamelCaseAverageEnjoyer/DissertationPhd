from numpy.random import uniform
from datetime import datetime
# import kiam_astro
from tiny_functions import * 
from plot_func import *
from cosmetic import *

SHAMANISM = {"KalmanQuaternionNormalize": True,   # Нормировка кватернионов в фильтре Калмана
             "KalmanSpinLimit": [True, 1e-3],  # Ограничение скорости вращения
             "ClohessyWiltshireC1=0": True}  # Траектории без дрейфа (без учёта аэродинамики)
GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas']
NAVIGATIONS = ['perfect', 'near', 'random']
R_V_CubeSat_SPREAD = [0, 0]
DISTORTION = 0.


def get_dipoles(r, ind: str):
    if len(ind) != 1 and ind not in "xyz":
        raise ValueError(f'Индекс "{ind}" должен быть 1 символом из [x, y, z]')
    cos_a = (r[0] * int(ind == 'x') + r[1] * int(ind == 'y') + r[2] * int(ind == 'z'))
    sin_a = (np.sqrt(r[1] ** 2 + r[2] ** 2) * int(ind == 'x') +
             np.sqrt(r[0] ** 2 + r[2] ** 2) * int(ind == 'y') +
             np.sqrt(r[0] ** 2 + r[1] ** 2) * int(ind == 'z'))
    aside = ((r[1] + r[2]) * int(ind == 'x') +
             (r[0] + r[2]) * int(ind == 'y') +
             r[2] * int(ind == 'z'))
    if abs(sin_a) > 1e-4:
        return np.cos(cos_a * np.pi / 2) / sin_a + DISTORTION * cos_a ** 2 + DISTORTION * aside
    return 1e-4 + DISTORTION * cos_a ** 2 + DISTORTION * aside

def get_gain(o, r: Union[float, np.ndarray], mode3: bool = False):
    r1 = r / np.linalg.norm(r)
    if o.gain_mode == GAIN_MODES[1]:
        return np.linalg.norm([r1[0] * 1, r1[1] * o.ellipsoidal_signal, r1[2] * o.ellipsoidal_signal])
    if o.gain_mode == GAIN_MODES[2]:
        return get_dipoles(r1, 'z')
    if o.gain_mode == GAIN_MODES[3] or (mode3 and o.gain_mode == GAIN_MODES[4]):
        return get_dipoles(r1, 'x') + get_dipoles(r1, 'y')
    if o.gain_mode == GAIN_MODES[4] and not mode3:
        return [get_dipoles(r1, 'x'), get_dipoles(r1, 'y')]
    return 1

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
        self.ellipsoidal_signal = 0.9
        self.gain_mode = GAIN_MODES[0]

        # Индивидуальные параметры
        self.r_orf = [uniform(-r_spread, r_spread, 3) for _ in range(self.n)]
        self.v_orf = [uniform(-v_spread, v_spread, 3) for _ in range(self.n)]
        self.w_orf = [uniform(-w_spread, w_spread, 3) for _ in range(self.n)]
        self.q, self.q_ = [[uniform(-1, 1, 4) for _ in range(self.n)] for _ in range(2)]
        self.c_hkw, self.line, self.line_kalman, self.line_difference, self.attitude_difference, self.spin_difference, \
            self.z_difference = [[[] for _ in range(self.n)] for _ in range(7)]
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_ = \
            [[[[] for _ in range(self.n)] for _ in range(self.n)] for _ in range(4)]
        for i in range(self.n):
            if SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
            self.q[i] /= np.linalg.norm(self.q[i])
            self.q_[i] /= np.linalg.norm(self.q_[i])
        prm_poor = [np.append(np.append(np.append(uniform(-r_spread, r_spread, 3), self.q_[i]),
                                        uniform(-v_spread, v_spread, 3)), uniform(-w_spread, w_spread, 3))
                    for i in range(self.n)]
        prm_good = [np.append(np.append(np.append(self.r_orf[i], self.q[i]), self.v_orf[i]), self.w_orf[i])
                    for i in range(self.n)]

        # Параметры начального приближения
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
        self.r_mass_center = np.array([uniform(-i, i) for i in self.mass_center_error])
        self.size = sizes[self.model_number]
        self.ellipsoidal_signal = 0.9
        self.gain_mode = GAIN_MODES[0]

        # Индивидуальные параметры
        self.r_orf = [uniform(-r_spread, r_spread, 3) for _ in range(self.n)]
        self.v_orf = [uniform(-v_spread, v_spread, 3) for _ in range(self.n)]
        self.w_orf = [uniform(-w_spread, w_spread, 3) for _ in range(self.n)]
        self.q = [np.array([uniform(-1, 1) for _ in range(4)]) for _ in range(self.n)]
        self.c_hkw, self.line = [[[] for _ in range(self.n)] for _ in range(2)]
        self.signal_power, self.real_dist, self.calc_dist, self.calc_dist_, self.kalm_dist = \
            [[[[] for _ in range(n_f)] for _ in range(self.n)] for _ in range(5)]
        for i in range(self.n):
            if SHAMANISM["ClohessyWiltshireC1=0"]:
                self.v_orf[i][0] = - 2 * self.r_orf[i][2] * w_orb
            self.q[i] /= np.linalg.norm(self.q[i])

        # Прорисовка ножек
        self.legs_x = 0.0085
        self.legs_z = 0.007

class KalmanFilter:
    def __init__(self, f: FemtoSat, c: CubeSat, p, dt: float, w_orb: float, kalman_coef: dict,
                 orientation: bool = False, single_femto_filter: bool = True):

        # Общие параметры
        self.f = f  # Фемтоспутники
        self.c = c  # Кубсаты
        self.p = p  # Динамическая модель
        self.dt = dt
        self.gain_mode = GAIN_MODES[0]
        self.q_coef = kalman_coef['q']
        self.p_coef = kalman_coef['p']
        self.r_matrix = kalman_coef['r']
        self.orientation = orientation
        self.single_femto_filter = single_femto_filter
        self.t = len(f.rv_orf_calc[0]) if orientation else 6
        self.sequence_under_diagonal = flatten([[[i, j] for i in range(j)] for j in range(self.f.n)])
        self.sigmas, self.real_sigmas = [[[] for _ in range(self.t * self.f.n)] for _ in range(2)]

        # Матрицы фильтра в начальный момент времени
        if not self.orientation:
            self.r_orf_estimation = np.array([np.append(f.rv_orf_calc[i][0:3],
                                                        f.rv_orf_calc[i][7:10]) for i in range(f.n)])
            self.phi_ = np.array([[1, 0, 0, dt, 0, 0],
                                  [0, 1, 0, 0, dt, 0],
                                  [0, 0, 1, 0, 0, dt],
                                  [0, 0, 0, 1, 0, -2 * w_orb * dt],
                                  [0, - w_orb ** 2 * dt, 0, 0, 1, 0],
                                  [0, 0, 3 * w_orb ** 2 * dt,  2 * w_orb * dt, 0, 1]])
            self.d_ = np.vstack([np.zeros((3, 3)), np.eye(3)])
            self.p_ = [np.diag([self.p_coef[0]]*3 + [self.p_coef[1]]*3) for _ in range(f.n)]
            self.q_ = np.diag([self.q_coef[0]]*3)
        else:  # 'rvq'
            self.r_orf_estimation = np.array(f.rv_orf_calc)
            self.phi_ = np.array([[1, 0, 0, 0, 0, 0, 0, dt, 0, 0, 0, 0, 0],
                                  [0, 1, 0, 0, 0, 0, 0, 0, dt, 0, 0, 0, 0],
                                  [0, 0, 1, 0, 0, 0, 0, 0, 0, dt, 0, 0, 0],
                                  [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 1, 0, -2 * w_orb * dt, 0, 0, 0],
                                  [0, -w_orb ** 2 * dt, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                  [0, 0, 3 * w_orb ** 2 * dt, 0, 0, 0, 0, 2 * w_orb * dt, 0, 1, 0, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]])
            self.p_ = [np.diag([self.p_coef[0]]*3 + [self.p_coef[2]]*4 + [self.p_coef[1]]*3 + [self.p_coef[3]]*3)
                       for _ in range(f.n)]
            # self.d_ = np.vstack([np.zeros((10, 6)), np.hstack([np.eye(3), np.eye(3)])])
            self.d_ = np.vstack([np.zeros((7, 6)),
                                 np.hstack([np.eye(3), np.zeros((3, 3))]),
                                 np.hstack([np.eye(3), np.eye(3)])])
            self.q_ = np.diag([self.q_coef[0]]*3 + [self.q_coef[1]]*3)

        # Расширешние на учёт несколько аппаратов в фильтре
        if not self.single_femto_filter:
            self.phi_ = np.bmat([[np.zeros([self.t, self.t])] * i + [self.phi_] +
                                 [np.zeros([self.t, self.t])] * (self.f.n - i - 1) for i in range(self.f.n)])
            tmp = 6 if self.orientation else 3
            self.q_ = np.bmat([[np.zeros([tmp, tmp])] * i + [self.q_] +
                               [np.zeros([tmp, tmp])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.p_ = np.bmat([[np.zeros([self.t, self.t])] * i + [self.p_[0]] +
                               [np.zeros([self.t, self.t])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.d_ = np.bmat([[np.zeros([self.t, tmp])] * i + [self.d_] +
                               [np.zeros([self.t, tmp])] * (self.f.n - i - 1) for i in range(self.f.n)])
            # self.d_ = np.bmat([[self.d_] for i in range(self.f.n)])

    def calc(self, i: int, i_c: int = 0) -> None:
        """
        :param i: ID-номер дочернего FemtoSat
        :param i_c: ID-номер дочернего CubeSat
        :return: None
        """
        z_ = self.c.calc_dist[0][i][-1]
        r_ = self.r_orf_estimation[i]
        if self.p.is_aero:
            rv_m = self.p.aero_integrate(self.f, i=i, r=r_[0:3], v=r_[7:10] if self.orientation else r_[3:6])
            if self.orientation:
                thw = self.p.attitude_integrate(self.f, i=i, th=r_[3:7], w=r_[10:13])
                r_m = np.append(np.append(np.append(rv_m[0], thw[0]), rv_m[1]), thw[1])
            else:
                r_m = np.append(rv_m[0], rv_m[1])
        else:
            r_m = self.phi_ @ r_
            if self.orientation:
                thw = self.p.attitude_integrate(self.f, i=i, th=r_[3:7], w=r_[10:13])
                r_m[3:7], r_m[10:13] = (thw[0], thw[1])

        if self.c.gain_mode == GAIN_MODES[4]:
            if self.orientation:
                tmp1 = get_gain(self.c, quart2dcm(self.c.q[i_c]) @ np.array(self.c.r_orf[i] - r_m[0:3]))
                tmp2 = get_gain(self.f, quart2dcm(r_m[3:7]) @ np.array(self.c.r_orf[i] - r_m[0:3]), mode3=True)
                signal_rate = [tmp1[ii] * tmp2 for ii in range(2)]
                z_model = np.array([np.linalg.norm(r_m[0:3] - self.c.r_orf[0]) / np.sqrt(signal_rate[i])
                                    for i in range(2)])
                z_ *= np.sqrt(signal_rate)
            else:
                z_model = np.array([1., 1.])
                z_ = np.append(z_, z_)

            self.f.z_difference[i] += [abs(z_model - z_)]

            h_ = np.array([[(r_m[j] - self.c.r_orf[0][j]) / z_model[ii] for j in range(3)] + (self.t - 3) * [0.]
                           for ii in range(2)])
            # print(f"Ф: {self.phi_.shape}, D: {self.d_.shape}, Q: {self.q_.shape}")
            q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt
            p_m = self.phi_ @ self.p_[i] @ self.phi_.T + q_tilda
            k_ = p_m @ h_.T @ np.linalg.pinv(h_ @ p_m @ h_.T + np.eye(2) * self.r_matrix)
            tmp = r_m + k_ @ (z_ - z_model)
            self.p_[i] = (np.eye(self.t) - k_ @ h_) @ p_m
        else:
            # print(f"r_m:{len(r_m)}, c.q:{len(self.c.q)}, r_orf:{len(self.c.r_orf)}, i={i}")
            signal_rate = get_gain(self.c, quart2dcm(self.c.q[i_c]) @ np.array(self.c.r_orf[i_c] - r_m[0:3])) * \
                get_gain(self.f, quart2dcm(r_m[3:7]) @ np.array(self.c.r_orf[i_c] - r_m[0:3])) \
                if self.orientation else 1
            z_model = np.linalg.norm(r_m[0:3] - self.c.r_orf[0])
            z_ *= np.sqrt(signal_rate)
            self.f.z_difference[i] += [abs(z_model - z_)]

            h_ = np.array([(r_m[j] - self.c.r_orf[0][j]) / z_model for j in range(3)] + (self.t - 3) * [0.])
            q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt
            p_m = self.phi_ @ self.p_[i] @ self.phi_.T + q_tilda
            k_ = p_m @ h_.T / (h_ @ p_m @ h_.T + self.r_matrix)
            tmp = r_m + k_ * (z_ - z_model)
            self.p_[i] = (np.eye(self.t) - np.outer(k_, h_)) @ p_m
        if SHAMANISM["KalmanQuaternionNormalize"] and self.orientation:
            tmp[3:7] = tmp[3:7] / np.linalg.norm(tmp[3:7])
        if SHAMANISM["KalmanSpinLimit"][0] and self.orientation and \
                np.linalg.norm(tmp[10:13]) > SHAMANISM["KalmanSpinLimit"][1]:
            tmp[10:13] = tmp[10:13] / np.linalg.norm(tmp[10:13]) * SHAMANISM["KalmanSpinLimit"][1]
        self.r_orf_estimation[i] = tmp

    def calc_all(self, i_c: int = 0) -> None:
        """
        :param i_c: ID-номер материнского CubeSat
        :return: None
        """
        ll = int(self.f.n * (self.f.n + 1) / 2)  # Рёбра полного графа N фемтосатов + 1 кубсат
        z_ = np.array([self.c.calc_dist[0][i][-1] for i in range(self.f.n)] + flatten([
                      [self.f.calc_dist[j][i][-1] for i in range(j)] for j in range(self.f.n)]))
        if self.p.is_aero:
            r_ = self.r_orf_estimation
            rv_m = [self.p.aero_integrate(self.f, i=i, r=r_[i][0:3], v=r_[i][7:10] if self.orientation else r_[i][3:6])
                    for i in range(self.f.n)]
            if self.orientation:
                thw = [self.p.attitude_integrate(self.f, i=i, th=r_[i][3:7], w=r_[i][10:13]) for i in range(self.f.n)]
                r_m = np.array(flatten([np.append(np.append(np.append(rv_m[i][0], thw[i][0]), rv_m[i][1]), thw[i][1])
                                        for i in range(self.f.n)]))
            else:
                r_m = np.array(flatten([np.append(rv_m[i][0], rv_m[i][1]) for i in range(self.f.n)]))
        else:
            r_ = np.array(flatten(self.r_orf_estimation))
            r_m = np.array(self.phi_ @ r_)[i_c]  # 6t_6t @ 6t -> 6t
            # r_m = self.phi_ @ r_
            if self.orientation:
                for i in range(self.f.n):
                    r__ = self.r_orf_estimation[i]
                    thw = self.p.attitude_integrate(self.f, i=i, th=r__[3:7], w=r__[10:13])
                    r_m[self.t*i+3:self.t*i+7], r_m[self.t*i+10:self.t*i+13] = (thw[0], thw[1])

        signal_rate = np.array([get_gain(self.f, quart2dcm(r_m[self.t*i+3:self.t*i+7]) @
                                         (r_m[self.t*i+0:self.t*i+3] - self.c.r_orf[i_c])) *
                                get_gain(self.c, quart2dcm(self.c.q[i_c]) @
                                         (r_m[self.t*i+0:self.t*i+3] - self.c.r_orf[i_c]))
                                for i in range(self.f.n)] +
                               flatten([[get_gain(self.f, quart2dcm(r_m[self.t*i+3:self.t*i+7])
                                                         @ (r_m[self.t*i+0:self.t*i+3] -
                                                            r_m[self.t*j+0:self.t*j+3])) *
                                         get_gain(self.f, quart2dcm(r_m[self.t*j+3:self.t*j+7])
                                                         @ (r_m[self.t*i+0:self.t*i+3] -
                                                            r_m[self.t*j+0:self.t*j+3]))
                                         for i in range(j)] for j in range(self.f.n)])) \
            if self.orientation else np.array([1] * int(self.f.n * (self.f.n + 1) // 2))

        z_model = np.array([np.linalg.norm(r_m[0+self.t*i:3+self.t*i] - self.c.r_orf[0]) for i in range(self.f.n)] +
                           flatten([[np.linalg.norm(r_m[0+self.t*j:3+self.t*j] - r_m[0+self.t*i:3+self.t*i])
                                     for i in range(j)] for j in range(self.f.n)]))
        z_ *= np.sqrt(signal_rate)
        for i in range(self.f.n):
            self.f.z_difference[i] += [abs(z_model[i] - z_[i])]

        def local_h_func(i: int, j: int) -> list:
            """Функция для построения матрицы H:
            i - строка из n + n(n-1)/2
            j - столбец из n (6 строк выдаёт как блок в return)"""
            r_k = r_m[self.t*j:self.t*j+3]
            if (i < self.f.n and i == j) or (i >= self.f.n and j in self.sequence_under_diagonal[i - self.f.n]):
                return [(r_k[j] - self.c.r_orf[0][j]) / z_model[i] for j in range(3)] + (self.t - 3) * [0.]
            else:
                return self.t * [0.]

        h_ = np.vstack([flatten([local_h_func(i, j) for i in range(self.f.n)]) for j in range(ll)])
        q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt  # nt_nt
        p_m = self.phi_ @ self.p_ @ self.phi_.T + q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        r = np.eye(ll) * self.r_matrix  # n + n(n-1)/2
        k_ = p_m @ h_.T @ np.linalg.inv(h_ @ p_m @ h_.T + r)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.p_ = (np.eye(self.t * self.f.n) - k_ @ h_) @ p_m
        r_orf_estimation = np.matrix(r_m) + k_ @ (z_ - z_model)
        for i in range(self.f.n):
            tmp = np.array(r_orf_estimation)[i_c][0+i*self.t:self.t+i*self.t]
            if SHAMANISM["KalmanQuaternionNormalize"]:
                tmp[3:7] = tmp[3:7] / np.linalg.norm(tmp[3:7])
            if SHAMANISM["KalmanSpinLimit"][0] and np.linalg.norm(tmp[10:13]) > SHAMANISM["KalmanSpinLimit"][1]:
                tmp[10:13] = tmp[10:13] / np.linalg.norm(tmp[10:13]) * SHAMANISM["KalmanSpinLimit"][1]
            self.r_orf_estimation[i] = tmp

class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat, kalman_coef: dict, method_navigation: str, dt: float = 1.,
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
        self.k = KalmanFilter(f=f, c=c, p=self, dt=self.dt, w_orb=self.w_orb, kalman_coef=kalman_coef,
                              orientation='rvq' in self.navigation, single_femto_filter='all' not in self.navigation)

        # Расчёт движения Хилла-Клохесси-Уилтшира
        self.c.c_hkw = [get_c_hkw(self.c.r_orf[i], self.c.v_orf[i], self.w_orb) for i in range(self.c.n)]
        self.f.c_hkw = [get_c_hkw(self.f.r_orf[i], self.f.v_orf[i], self.w_orb) for i in range(self.f.n)]

    # Интегрирование движения
    def aero_integrate(self, obj, i: int, r=None, v=None, th=None, w=None):
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
        return self.rk4_translate(r, v, self.get_full_acceleration(c_resist=c_resist, rho=rho, square=square,
                                                                   m=obj.mass, r=r, v=v))

    def attitude_integrate(self, obj, i: int, th=None, w=None):
        th = obj.q[i] if th is None else th
        w = obj.w_orf[i] if w is None else w
        return self.rk4_attitude(th, w, self.get_full_torgue(th=th, w=w))

    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.dt

        # Утверждение параметров
        if self.iter == 1:
            tmp1 = f", сферичность={int(self.c.ellipsoidal_signal*100)}%" if self.c.gain_mode == GAIN_MODES[1] \
                else ""
            tmp2 = f", сферичность={int(self.f.ellipsoidal_signal*100)}%" if self.f.gain_mode == GAIN_MODES[1] \
                else ""
            tmp = ", ориентации\n" if 'rvq' in self.navigation else "\n"
            my_print(f"Диаграмма антенн кубсата: {self.c.gain_mode}{tmp1}\n"
                     f"Диаграмма антенн фемтосатов: {self.f.gain_mode}{tmp2}\n"
                     f"Учёт аэродинамики: {self.is_aero}\n"
                     f"Применяется фильтр Калмана для поправки: положений, скоростей{tmp}" 
                     f"Фильтр Калмана основан на: "
                     f"{'одном чипсате' if self.k.single_femto_filter else 'всех чипсатах'}", color='c')

        # Вращательное движение
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                obj.q[i], obj.w_orf[i] = self.attitude_integrate(obj, i)

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
            for i_f in range(self.f.n):
                r_tmp = self.f.r_orf[i_f] - self.c.r_orf[i_c]
                self.c.real_dist[i_c][i_f] += [np.linalg.norm(r_tmp)]
                A_f = quart2dcm(np.array(self.f.q[i_f]))
                A_c = quart2dcm(np.array(self.c.q[i_c]))
                calc_dist = np.linalg.norm(r_tmp) / np.sqrt(get_gain(self.f, r=A_f @ r_tmp, mode3=True) *
                                                            get_gain(self.c, r=A_c @ r_tmp, mode3=True)) + \
                    np.random.normal(0, np.sqrt(self.r_matrix))
                self.c.calc_dist[i_c][i_f] += [calc_dist]
                self.c.kalm_dist[i_c][i_f] += [np.linalg.norm(self.f.r_orf[i_f] -
                                                              self.k.r_orf_estimation[i_f][0:3])]
                if self.iter % self.show_rate == 0:
                    self.f.line_kalman[i_f] += [self.k.r_orf_estimation[i_f][0],
                                                self.k.r_orf_estimation[i_f][1],
                                                self.k.r_orf_estimation[i_f][2]]
                    self.f.line_difference[i_f] += \
                        [np.array(self.k.r_orf_estimation[i_f][0:3] - np.array(self.f.r_orf[i_f]))]
                    if self.k.orientation:
                        self.f.attitude_difference[i_f] += [self.k.r_orf_estimation[i_f][3:7]
                                                            - np.array(self.f.q[i_f])]
                        self.f.spin_difference[i_f] += [self.k.r_orf_estimation[i_f][10:13]
                                                        - np.array(self.f.w_orf[i_f])]
        for i_f1 in range(self.f.n):
            for i_f2 in range(self.f.n):
                if i_f1 != i_f2:
                    r_tmp = self.f.r_orf[i_f1] - self.f.r_orf[i_f2]
                    self.f.real_dist[i_f1][i_f2] += [np.linalg.norm(r_tmp)]

                    q1 = self.f.q[i_f1]
                    q2 = self.f.q[i_f1]
                    calc_dist = np.linalg.norm(r_tmp) / np.sqrt(get_gain(self.f, r=quart2dcm(q1) @ r_tmp, mode3=True) *
                                                                get_gain(self.f, r=quart2dcm(q2) @ r_tmp, mode3=True)) \
                                + np.random.normal(0, np.sqrt(self.r_matrix))
                    self.f.calc_dist[i_f1][i_f2] += [calc_dist]
                else:
                    self.f.calc_dist[i_f1][i_f2] += [0.]

        # Оценка положения фемтоспутников
        if self.k.single_femto_filter:
            for i in range(self.f.n):
                self.k.calc(i)
                if self.c.gain_mode != GAIN_MODES[4]:
                    for j_n in range(self.f.n):
                        for j_t in range(9 if self.k.orientation else 3):  # range(self.k.t)
                            tmp = np.append(np.append(self.f.r_orf[j_n], self.f.v_orf[j_n]), list(self.f.q[j_n][1:4]))
                            self.k.sigmas[j_n * self.k.t + j_t] += [np.sqrt(self.k.p_[j_n][j_t][j_t]) * tmp[j_t]]
        else:
            self.k.calc_all()

    def integrate(self, t: float, animate: bool = False) -> None:
        my_print(f"Оборотов вокруг Земли: {round(10 * t / (3600 * 1.5)) / 10}  "
                 f"(дней: {round(100 * t / (3600 * 24)) / 100})", color='b')
        n = int(t // self.dt)
        flag = [0., 0.]
        for i in range(n):
            if i / n > (flag[0] + 0.1):
                flag[0] += 0.1
                per = int(10 * i / n)
                my_print(f"{10 * per}% [{'#' * per + ' ' * (10 - per)}]" +
                         real_workload_time(n=per, n_total=10, time_begin=self.time_begin,
                                            time_now=datetime.now()), color='m')
            if animate and i / n > (flag[1] + 0.01):
                flag[1] += 0.01
                plot_all(self, save=True, count=int(flag[1] // 0.01))
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

    def get_full_torgue(self, th, w) -> np.ndarray:
        return np.zeros(3)

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

    def rk4_attitude(self, L, w, e):
        def lw_right_part(Lw1, e1):
            L1, w1 = Lw1[0:4], Lw1[4:7]
            dL = 1 / 2 * q_dot([0, w1[0], w1[1], w1[2]], L1)
            return np.append(dL, e1)
        a = len(L)
        L_ = L if a == 4 else vec2quat(L)
        Lw = np.append(L_, w)
        k1 = lw_right_part(Lw, e)
        k2 = lw_right_part(Lw + k1 * self.dt / 2, e)
        k3 = lw_right_part(Lw + k2 * self.dt / 2, e)
        k4 = lw_right_part(Lw + k3 * self.dt, e)
        Lw = self.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        tmp = (Lw[0:4] + L_)  # / np.linalg.norm(Lw[0:4] + L_)
        return tmp[4-a:4], Lw[a:a+3] + w

class Objects:
    def __init__(self, method_navigation: str = None, kalman_coef: dict = None, dt: float = 1,
                 n_c: int = 1, n_f: int = 5, model_c: str = '1U',
                 start_navigation: str = 'perfect', start_navigation_tolerance: float = 0.9):
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
        self.c = CubeSat(n=n_c, n_f=n_f, model=model_c, w_orb=self.w_orb)
        self.f = FemtoSat(n=n_f, start_navigation_tolerance=start_navigation_tolerance,
                          start_navigation=start_navigation, w_orb=self.w_orb)
        self.p = PhysicModel(c=self.c, f=self.f, dt=dt,
                             method_navigation='kalman_filter rv' if method_navigation is None else method_navigation,
                             kalman_coef={'q': 1e-12, 'p': [1e-8] * 3, 'r': 1e-1}
                             if kalman_coef is None else kalman_coef, h_orb=self.h_orb)
