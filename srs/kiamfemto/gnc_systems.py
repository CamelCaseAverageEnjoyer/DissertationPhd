from spacecrafts import *
from data.observability_mapping_partial_derivatives import *

# >>>>>>>>>>>> Guidance <<<<<<<<<<<<
def guidance(c: CubeSat, f: FemtoSat, v: Variables, earth_turn: float) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            pass
            # if obj.operating_modes[i] == v.OPERATING_MODES_CHANGE[1]:  # Отсутствие аккумулятора на чипсате
            #     if earth_turn % 1 < 0.5 and obj.operating_mode[i] == v.OPERATING_MODES[-1]:
            #         obj.operating_mode[i] = v.OPERATING_MODES[0]
            #     if earth_turn % 1 > 0.5 and obj.operating_mode[i] != v.OPERATING_MODES[-1]:
            #         obj.operating_mode[i] = v.OPERATING_MODES[-1]

# >>>>>>>>>>>> Navigation <<<<<<<<<<<<
class KalmanFilter:
    """Оцениваемые параметры: ['r orf', 'q-3 irf', 'v orf', 'w irf']; согласовано с spacecrafts.py"""
    def __init__(self, f: FemtoSat, c: CubeSat, p: any):
        # Общие параметры
        self.f = f  # Фемтоспутники
        self.c = c  # Кубсаты
        self.p = p  # Динамическая модель
        self.v = p.v
        # self.sequence_under_diagonal = flatten([[[i, j] for i in range(j)] for j in range(self.f.n)])  # Удалить?

        self.estimation_params = self.params_dict2vec(d=f.apriori_params, separate_spacecraft=True)
        self.j = len(self.estimation_params[0])  # Вектор состояния 1 чипсата

        # Матрицы фильтра в начальный момент времени
        w0 = self.v.W_ORB
        if not self.v.NAVIGATION_ANGLES:  # Вектор состояния содержит только положение и скорость
            self.Phi = np.array([[0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 1, 0],
                                 [0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, -2*w0],
                                 [0, - w0**2, 0, 0, 0, 0],
                                 [0, 0, 3*w0**2,  2*w0, 0, 0]]) * self.v.dT + np.eye(self.j)
            self.D = np.vstack([np.zeros((3, 3)), np.eye(3)])
            self.P = [np.diag([self.v.KALMAN_COEF['p'][0]] * 3 + [self.v.KALMAN_COEF['p'][1]] * 3) for _ in range(f.n)]
            self.Q = np.diag([self.v.KALMAN_COEF['q'][0]] * 3)
        else:  # Вектор состояния содержит ещё и угловые переменные
            self.Phi = np.array([[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1/2],
                                 [0, 0, 0, 0, 0, 0, 0, 0, -2*w0, 0, 0, 0],
                                 [0, -w0**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 3*w0**2, 0, 0, 0, 2*w0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) * self.v.dT + np.eye(self.j)
            self.P = [np.diag([self.v.KALMAN_COEF['p'][0]] * 3 + [self.v.KALMAN_COEF['p'][2]] * 3 +
                              [self.v.KALMAN_COEF['p'][1]] * 3 + [self.v.KALMAN_COEF['p'][3]] * 3) for _ in range(f.n)]
            self.D = np.vstack([np.zeros((6, 6)), np.eye(6)])
            self.Q = np.diag([self.v.KALMAN_COEF['q'][0]] * 3 + [self.v.KALMAN_COEF['q'][1]] * 3)

        # Расширешние на учёт несколько аппаратов в фильтре
        if self.v.NAVIGATION_BY_ALL:
            # Тут написан очень интересный способ сделать блочно-диагональные-провальные-нахальные матрицы
            tmp = 6 if self.v.NAVIGATION_ANGLES else 3
            self.Phi = np.bmat([[np.zeros([self.j, self.j])] * i + [self.Phi] +
                                [np.zeros([self.j, self.j])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.Q = np.bmat([[np.zeros([tmp, tmp])] * i + [self.Q] +
                              [np.zeros([tmp, tmp])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.P = np.bmat([[np.zeros([self.j, self.j])] * i + [self.P[i]] +
                              [np.zeros([self.j, self.j])] * (self.f.n - i - 1) for i in range(self.f.n)])
            self.D = np.bmat([[np.zeros([self.j, tmp])] * i + [self.D] +
                              [np.zeros([self.j, tmp])] * (self.f.n - i - 1)
                              for i in range(self.f.n)])

        # Вывод
        if self.v.IF_ANY_PRINT:
            my_print(f"Матрицы Ф:{self.Phi.shape}, Q:{self.Q.shape}, P:{self.P.shape}, D:{self.D.shape}", color='g')

    def params_vec2dict(self, params: list = None, j: int = None, separate_spacecraft: bool = True):
        p = self.estimation_params if params is None else params
        j = self.j if j is None else j
        if self.v.NAVIGATION_ANGLES:
            if separate_spacecraft:
                r_orf = [p[i][0: 3] for i in range(self.f.n)]
                q_irf = [p[i][3: 6] for i in range(self.f.n)]
                v_orf = [p[i][6: 9] for i in range(self.f.n)]
                w_irf = [p[i][9: 12] for i in range(self.f.n)]
            else:
                r_orf = [p[i*j + 0: i*j + 3] for i in range(self.f.n)]
                q_irf = [p[i*j + 3: i*j + 6] for i in range(self.f.n)]
                v_orf = [p[i*j + 6: i*j + 9] for i in range(self.f.n)]
                w_irf = [p[i*j + 9: i*j + 12] for i in range(self.f.n)]
        else:
            if separate_spacecraft:
                r_orf = [p[i][0: 3] for i in range(self.f.n)]
                v_orf = [p[i][3: 6] for i in range(self.f.n)]
            else:
                r_orf = [p[i*j + 0: i*j + 3] for i in range(self.f.n)]
                v_orf = [p[i*j + 3: i*j + 6] for i in range(self.f.n)]
            q_irf, w_irf = [[None for _ in range(self.f.n)] for _ in range(2)]
        return {'r orf': r_orf, 'v orf': v_orf, 'w irf': w_irf, 'q-3 irf': q_irf}

    def params_dict2vec(self, d: dict, separate_spacecraft: bool = True):
        variables = ['r orf', 'q-3 irf', 'v orf', 'w irf'] if self.v.NAVIGATION_ANGLES else ['r orf', 'v orf']
        if separate_spacecraft:  # Вроде как, нужен только этот вариант
            return [np.array([d[v][i][j] for v in variables for j in range(3)]) for i in range(self.f.n)]
        else:
            return np.array([d[v][i][j] for i in range(self.f.n) for v in variables for j in range(3)])

    def get_estimation(self, i_f: int, v: str):
        d = self.params_vec2dict()
        return d[v][i_f]

    def calc(self) -> None:  # Считается, что NAVIGATION_BY_ALL = True
        from primary_info import measure_antennas_power
        from dynamics import rk4_translate, rk4_attitude

        # >>>>>>>>>>>> Предварительный расчёт <<<<<<<<<<<<
        c_take_len = len(get_gain(v=self.v, obj=self.c, r=np.ones(3), if_take=True))
        c_send_len = len(get_gain(v=self.v, obj=self.c, r=np.ones(3), if_send=True))
        f_take_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_take=True))
        f_send_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_send=True))
        z_len = int(self.f.n * self.c.n * (c_take_len*f_send_len + c_send_len*f_take_len) +
                    self.f.n * (self.f.n - 1) * f_take_len*f_send_len)
        my_print(f"Вектор измерений z_len: {z_len}     =     {self.f.n}*{self.c.n}*({c_take_len}*{f_send_len}+"
                 f"{c_send_len}*{f_take_len}) + {self.f.n}({self.f.n}-1){f_take_len}*{f_send_len}",
                 color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)

        # >>>>>>>>>>>> Этап экстраполяции <<<<<<<<<<<<
        d = self.params_vec2dict()

        # Моделирование орбитального движения на dT -> вектор состояния x_m
        if np.any(self.v.DYNAMIC_MODEL.values()):  # Если J2 или aero drag
            rv_m = [rk4_translate(v_=self.v, obj=self.f, i=i, r=d['r orf'][i], v=d['v orf'][i])
                    for i in range(self.f.n)]
            x_m = self.params_dict2vec(d={'r orf': [rv_m[i][0] for i in range(self.f.n)],
                                          'v orf': [rv_m[i][1] for i in range(self.f.n)],
                                          'q-3 irf': np.zeros(3), 'w irf': np.zeros(3)}, separate_spacecraft=False)
            # tmp = int(self.v.NAVIGATION_ANGLES)
            # x_m = np.array(flatten([np.append(np.append(np.append(rv_m[i][0], [0]*4*tmp), rv_m[i][1]), [0]*3*tmp)
            #                         for i in range(self.f.n)]))
        else:
            x_m = np.array(self.Phi @ np.array(flatten(self.estimation_params)))[0]

        # Моделирование углового движения на dT -> вектор состояния x_m
        if self.v.NAVIGATION_ANGLES:
            tmp = self.params_vec2dict(x_m)
            for i in range(self.f.n):
                tmp['q-3 irf'][i], tmp['w irf'][i] = rk4_attitude(v_=self.v, obj=self.f, i=i,
                                                                  q=d['q-3 irf'][i], w=d['w irf'][i])
            x_m = self.params_dict2vec(tmp)

        # Измерения с поправкой на угловой коэффициент усиления G (signal_rate)
        dr_cf, dr_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v, get_ready_measurements=True)
        my_print(f"Длина длин 1: {len(dr_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        my_print(f"Длина длин 2: {len(dr_cf + dr_ff)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        z_ = np.array(dr_cf + dr_ff)
        # Поправка на G ниже
        if self.v.NAVIGATION_ANGLES:
            signal_rate_cf, signal_rate_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v,
                                                                           get_signal_rates=True, j=self.j, x_m=x_m)
            my_print(f"Длина G 1: {len(signal_rate_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
            my_print(f"Длина G 2: {len(signal_rate_cf + signal_rate_ff)}", color='r',
                     if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
            signal_rate = np.array(signal_rate_cf + signal_rate_ff)
        else:
            signal_rate = np.array([1] * z_len)
        z_ *= np.sqrt(signal_rate)

        # Измерения согласно модели
        z_cf, z_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v, get_model_state=True, j=self.j,
                                                   x_m=x_m)
        my_print(f"Длина модельных длин 1: {len(z_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        my_print(f"Длина модельных длин 2: {len(z_cf + z_ff)}", color='r',
                 if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        z_model = np.array([z_cf + z_ff])

        for i_c in range(self.c.n):  # ПРОВЕРИТЬ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            for i_f in range(self.f.n):
                tmp = c_take_len * f_send_len + c_send_len * f_take_len
                # print(f"{i_c}:{i_f}   1 ---> {z_cf[i_f*tmp: (i_f+1)*tmp]}")
                # print(f"{i_c}:{i_f}   2 ---> {z_[i_f*tmp: (i_f+1)*tmp]}")
                self.c.z_difference[i_c][i_f] += [np.linalg.norm(np.array(z_cf[i_f*tmp: (i_f+1)*tmp]) -
                                                                 np.array(z_[i_f*tmp: (i_f+1)*tmp]))]

        # Расчёт матриц
        Q_tilda = self.Phi @ self.D @ self.Q @ self.D.T @ self.Phi.T  # * self.v.dT  # nt_nt
        P_m = self.Phi @ self.P @ self.Phi.T + Q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        # >>>>>>>>>>>> Этап коррекции <<<<<<<<<<<<
        # ПОМЕНЯТЬ f.q, c.q НА ПРЕДПОЛАГАЕМЫЕ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        q_f = [self.f.q[i][1:4] for i in range(self.f.n)]
        q_c = [self.c.q[i][1:4] for i in range(self.c.n)]
        H, _, _ = h_matrix(c_ant=self.v.N_ANTENNA_C, f_ant=self.v.N_ANTENNA_F, fn=self.f.n, cn=self.c.n,
                           angles_navigation=self.v.NAVIGATION_ANGLES, r_f=self.f.r_orf, r_c=self.c.r_orf,
                           multy_antenna_send=self.v.MULTI_ANTENNA_SEND, multy_antenna_take=self.v.MULTI_ANTENNA_TAKE,
                           w_0=self.v.W_ORB, t=self.p.t, q_f=q_f, q_c=q_c)

        R = np.eye(z_len) * self.v.KALMAN_COEF['r']  # n + n(n-1)/2
        k_ = P_m @ H.T @ np.linalg.inv(H @ P_m @ H.T + R)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.P = (np.eye(self.j * self.f.n) - k_ @ H) @ P_m
        r_orf_estimation = np.array((np.matrix(x_m).T + k_ @ (z_ - z_model).T).T)[0]
        for i in range(self.f.n):
            tmp = r_orf_estimation[(0 + i) * self.j: (1 + i) * self.j]
            if self.v.SHAMANISM["KalmanQuaternionNormalize"] and self.v.NAVIGATION_ANGLES:
                tmp2 = vec2quat(tmp[3:6])
                tmp[3:6] = (tmp2 / np.linalg.norm(tmp2))[1:4]
            if self.v.SHAMANISM["KalmanSpinLimit"][0] and \
                    np.linalg.norm(tmp[9:12]) > self.v.SHAMANISM["KalmanSpinLimit"][1]:
                tmp[9:12] = tmp[9:12] / np.linalg.norm(tmp[9:12]) * self.v.SHAMANISM["KalmanSpinLimit"][1]
            if self.v.SHAMANISM["KalmanVelocityLimit"][0] and \
                    np.linalg.norm(tmp[6:9]) > self.v.SHAMANISM["KalmanVelocityLimit"][1]:
                tmp[6:9] = tmp[6:9] / np.linalg.norm(tmp[6:9]) * self.v.SHAMANISM["KalmanSpinLimit"][1]
            if self.v.SHAMANISM["KalmanPositionLimit"][0] and \
                    np.linalg.norm(tmp[0:3]) > self.v.SHAMANISM["KalmanPositionLimit"][1]:
                tmp[0:3] = tmp[0:3] / np.linalg.norm(tmp[0:3]) * self.v.SHAMANISM["KalmanPositionLimit"][1]
            self.estimation_params[i] = tmp

        # Запись и отображение
        if self.p.iter == 1:
            with open("kiamfemto/data/measures_vector_notes_last.txt", "w") as f:
                f.write("# Рассчитано в PyCharm\n# Параметры: {rel} {N_1} {N_2} {i_1} {i_2} {send_len} {take_len}\n")
                f.write(f"# Параметр CUBESAT_AMOUNT {self.v.CUBESAT_AMOUNT}\n")
                f.write(f"# Параметр CHIPSAT_AMOUNT {self.v.CHIPSAT_AMOUNT}\n")
                f.write(f"# Параметр NAVIGATION_ANGLES {self.v.NAVIGATION_ANGLES}\n")
                f.write(f"# Параметр N_ANTENNA_C {self.v.N_ANTENNA_C}\n")
                f.write(f"# Параметр N_ANTENNA_F {self.v.N_ANTENNA_F}\n")
                f.write(f"# Параметр MULTI_ANTENNA_TAKE {self.v.MULTI_ANTENNA_TAKE}\n")
                f.write(f"# Параметр MULTI_ANTENNA_SEND {self.v.MULTI_ANTENNA_SEND}\n")
                for j in range(z_len):
                    f.write(f"{self.v.MEASURES_VECTOR_NOTES[j]}\n")
            my_print(f"P-: {P_m.shape}, H.T: {H.T.shape}, H: {H.shape}, R: {R.shape}", color='g')
            my_print(f"x-: {np.matrix(x_m).shape}, K: {k_.shape}, z: {(z_ - z_model).shape}", color='g')
            my_print(f"r_orf_estimation: {r_orf_estimation.shape}", color='g')
            # my_print(H, color='y')
            # my_print(self.Phi, color='y')

def navigate(k: KalmanFilter):
    k.calc()  # Пока что при любом OPERATING_MODES (если весь рой выпал)

# >>>>>>>>>>>> Control <<<<<<<<<<<<
