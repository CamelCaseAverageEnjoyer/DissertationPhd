from srs.kiamfemtosat.spacecrafts import *

# >>>>>>>>>>>> Guidance <<<<<<<<<<<<
def guidance(c: CubeSat, f: FemtoSat, earth_turn: float) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            if obj.operating_modes[i] == OPERATING_MODES_CHANGE[1]:  # Отсутствие аккумулятора на чипсате
                if earth_turn % 1 < 0.5 and obj.operating_mode[i] == OPERATING_MODES[-1]:
                    obj.operating_mode[i] = OPERATING_MODES[0]
                if earth_turn % 1 > 0.5 and obj.operating_mode[i] != OPERATING_MODES[-1]:
                    obj.operating_mode[i] = OPERATING_MODES[-1]

# >>>>>>>>>>>> Navigation <<<<<<<<<<<<
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
        self.t = len(f.rv_orf_calc[0]) if orientation else 6  # Вектор состояния 1 чипсата?
        self.sequence_under_diagonal = flatten([[[i, j] for i in range(j)] for j in range(self.f.n)])
        self.sigmas, self.real_sigmas = [[[] for _ in range(self.t * self.f.n)] for _ in range(2)]

        # Матрицы фильтра в начальный момент времени
        if not NAVIGATION_ANGLES:
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
        if NAVIGATION_BY_ALL:
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

    def new_calc(self) -> None:  # Считается, что NAVIGATION_BY_ALL = True
        z_len = int(self.f.n * self.c.n) + int(self.f.n * (self.f.n - 1) / 2)  # ll = Рёбра НЕполного графа
        # Расчёт. Примечание: вектор _self.r_orf_estimation_ должен быть длины _self.f.n_ * _self.t_ (вроде)
        r_ = self.r_orf_estimation  # [i_f][i_t]
        if AERO_DRAG:  # Переделать весь расчёт. Это какой-то крокодил и он тебя сожрёт рано или поздно
            rv_m = [
                self.p.aero_integrate(self.f, i=i, r=r_[i][0:3], v=r_[i][7:10] if self.orientation else r_[i][3:6])
                for i in range(self.f.n)]
            if NAVIGATION_ANGLES:
                thw = [self.p.attitude_integrate(self.f, i=i, th=r_[i][3:7], w=r_[i][10:13]) for i in
                       range(self.f.n)]
                r_m = np.array(
                    flatten([np.append(np.append(np.append(rv_m[i][0], thw[i][0]), rv_m[i][1]), thw[i][1])
                             for i in range(self.f.n)]))
            else:
                r_m = np.array(flatten([np.append(rv_m[i][0], rv_m[i][1]) for i in range(self.f.n)]))
        else:
            r_m = np.array(self.phi_ @ np.array(flatten(r_)))[0]  # self.f.n * self.t | flatten -> [1..t] + [1..t] ...
            if NAVIGATION_ANGLES:
                for i in range(self.f.n):
                    r__ = self.r_orf_estimation[i]
                    q_w = self.p.attitude_integrate(self.f, i=i, th=r__[3:7], w=r__[10:13])
                    r_m[i*self.t + 3:i*self.t + 7], r_m[i*self.t + 10:i*self.t + 13] = (q_w[0], q_w[1])

        # Измерения с поправкой на signal_rate
        z_ = np.array([])
        for i_c in range(self.c.n):
            z_ = np.append(z_, [self.c.calc_dist[i_c][i][-1] for i in range(self.f.n)])
        z_ = np.append(z_, flatten([[self.f.calc_dist[j][i][-1] for i in range(j)] for j in range(self.f.n)]))

        signal_rate = [] if NAVIGATION_ANGLES else np.array([1] * z_len)
        if NAVIGATION_ANGLES:
            for i_c in range(self.c.n):
                signal_rate = np.append(signal_rate,
                                        [get_gain(self.f, quart2dcm(r_m[i*self.t + 3:i*self.t + 7]) @
                                                  (r_m[i*self.t + 0:i*self.t + 3] - self.c.r_orf[i_c]))[0] *
                                         get_gain(self.c, quart2dcm(self.c.q[i_c]) @
                                                  (r_m[i*self.t + 0:i*self.t + 3] - self.c.r_orf[i_c]))[0]
                                         for i in range(self.f.n)])
            signal_rate = np.append(signal_rate, flatten(
                [[get_gain(self.f, quart2dcm(r_m[i*self.t + 3:i*self.t + 7]) @
                           (r_m[i*self.t + 0:i*self.t + 3] - r_m[j*self.t + 0:j*self.t + 3]))[0] *
                  get_gain(self.f, quart2dcm(r_m[j*self.t + 3:j*self.t + 7]) @
                           (r_m[i*self.t + 0:i*self.t + 3] - r_m[j*self.t + 0:j*self.t + 3]))[0]
                  for i in range(j)] for j in range(self.f.n)]))
        z_ *= np.sqrt(signal_rate)

        # Измерения согласно модели
        z_model = np.array([])
        for i_c in range(self.c.n):
            z_model = np.append(z_model, [np.linalg.norm(r_m[i*self.t + 0:i*self.t + 3] - self.c.r_orf[0])
                                          for i in range(self.f.n)])
        z_model = np.append(z_model,
                            flatten([[np.linalg.norm(r_m[j*self.t + 0:j*self.t + 3] - r_m[i*self.t + 0:i*self.t + 3])
                                      for i in range(j)] for j in range(self.f.n)]))
        for i in range(self.f.n):
            self.f.z_difference[i] += [abs(z_model[i] - z_[i])]

        def local_h_func(i: int, j: int) -> list:
            """Функция для построения матрицы H:
            i - блок столбцов из f.n
            j - строка из n (self.t строк выдаёт как блок в return)
            сперва int(self.f.n * self.c.n) строк
            потом  int(self.f.n * (self.f.n - 1) / 2) строк"""
            r_k = r_m[i*self.t + 0:i*self.t + 3]
            if (j < self.c.n * self.f.n and i == j) or \
                    (i >= self.c.n * self.f.n and j in self.sequence_under_diagonal[i - self.c.n * self.f.n]):
                return [(r_k[jj] - self.c.r_orf[i][jj]) / z_model[i] for jj in range(3)] + (self.t - 3) * [0.]
            else:
                return self.t * [0.]

        h_ = np.vstack([flatten([local_h_func(i, j) for i in range(self.f.n)]) for j in range(z_len)])
        q_tilda = self.phi_ @ self.d_ @ self.q_ @ self.d_.T @ self.phi_.T * self.dt  # nt_nt
        p_m = self.phi_ @ self.p_ @ self.phi_.T + q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        r = np.eye(z_len) * self.r_matrix  # n + n(n-1)/2
        k_ = p_m @ h_.T @ np.linalg.inv(h_ @ p_m @ h_.T + r)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.p_ = (np.eye(self.t * self.f.n) - k_ @ h_) @ p_m
        r_orf_estimation = np.matrix(r_m) + k_ @ (z_ - z_model)
        for i in range(self.f.n):
            tmp = np.array(r_orf_estimation)[0 + i * self.t:self.t + i * self.t]
            if SHAMANISM["KalmanQuaternionNormalize"] and NAVIGATION_ANGLES:
                tmp[3:7] = tmp[3:7] / np.linalg.norm(tmp[3:7])
            if SHAMANISM["KalmanSpinLimit"][0] and np.linalg.norm(tmp[10:13]) > SHAMANISM["KalmanSpinLimit"][1]:
                tmp[10:13] = tmp[10:13] / np.linalg.norm(tmp[10:13]) * SHAMANISM["KalmanSpinLimit"][1]
            self.r_orf_estimation[i] = tmp

    def calc(self, i: int, i_c: int = 0) -> None:
        """
        :param i: ID-номер дочернего FemtoSat
        :param i_c: ID-номер материнского CubeSat
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
            signal_rate = get_gain(self.c, quart2dcm(self.c.q[i_c]) @ np.array(self.c.r_orf[i_c] - r_m[0:3]))[0] * \
                          get_gain(self.f, quart2dcm(r_m[3:7]) @ np.array(self.c.r_orf[i_c] - r_m[0:3]))[0] \
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
                                         (r_m[self.t*i+0:self.t*i+3] - self.c.r_orf[i_c]))[0] *
                                get_gain(self.c, quart2dcm(self.c.q[i_c]) @
                                         (r_m[self.t*i+0:self.t*i+3] - self.c.r_orf[i_c]))[0]
                                for i in range(self.f.n)] +
                               flatten([[get_gain(self.f, quart2dcm(r_m[self.t*i+3:self.t*i+7])
                                                  @ (r_m[self.t*i+0:self.t*i+3] - r_m[self.t*j+0:self.t*j+3]))[0] *
                                         get_gain(self.f, quart2dcm(r_m[self.t*j+3:self.t*j+7])
                                                  @ (r_m[self.t*i+0:self.t*i+3] - r_m[self.t*j+0:self.t*j+3]))[0]
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

# >>>>>>>>>>>> Control <<<<<<<<<<<<
