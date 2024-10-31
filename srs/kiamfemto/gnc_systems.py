import numpy as np

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
            i_f = 0
            w = f.apriori_params['w irf'][i_f]
            Phi_w = np.linalg.inv(f.J) @ (- get_antisymmetric_matrix(w) @ f.J + get_antisymmetric_matrix(f.J @ w))
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
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) * self.v.dT + np.eye(self.j) + \
                       np.bmat([[np.zeros((9, 12))], [np.zeros((3, 9)), Phi_w]])
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
        if True:
            f = self.f
            c = self.c
            v = self.v
            p = self.p
            j = self.j
            c_take_len = len(get_gain(v=v, obj=c, r=np.ones(3), if_take=True))
            c_send_len = len(get_gain(v=v, obj=c, r=np.ones(3), if_send=True))
            f_take_len = len(get_gain(v=v, obj=f, r=np.ones(3), if_take=True))
            f_send_len = len(get_gain(v=v, obj=f, r=np.ones(3), if_send=True))
            z_len = int(f.n * c.n * (c_take_len*f_send_len + c_send_len*f_take_len) +
                        f.n * (f.n - 1) * f_take_len*f_send_len)
            my_print(f"Вектор измерений z_len: {z_len}     =     {f.n}*{c.n}*({c_take_len}*{f_send_len}+"
                     f"{c_send_len}*{f_take_len}) + {f.n}({f.n}-1){f_take_len}*{f_send_len}",
                     color='r', if_print=p.iter == 1 and v.IF_TEST_PRINT)

        # >>>>>>>>>>>> Этап экстраполяции <<<<<<<<<<<<
        d = self.params_vec2dict()

        # Моделирование орбитального движения на dT -> вектор состояния x_m
        if np.any(v.DYNAMIC_MODEL.values()):  # Если J2 или aero drag
            rv_m = [rk4_translate(v_=v, obj=f, i=i, r=d['r orf'][i], v=d['v orf'][i])
                    for i in range(f.n)]
            x_m = self.params_dict2vec(d={'r orf': [rv_m[i][0] for i in range(f.n)],
                                          'v orf': [rv_m[i][1] for i in range(f.n)],
                                          'q-3 irf': np.zeros((f.n, 3)), 'w irf': np.zeros((f.n, 3))},
                                       separate_spacecraft=False)
        else:
            x_m = np.array(self.Phi @ np.array(flatten(self.estimation_params)))[0]

        # Моделирование углового движения на dT -> вектор состояния x_m
        if v.NAVIGATION_ANGLES:
            tmp = self.params_vec2dict(params=x_m, separate_spacecraft=False)
            for i in range(f.n):
                tmp['q-3 irf'][i], tmp['w irf'][i] = rk4_attitude(v_=v, obj=f, i=i,
                                                                  q=d['q-3 irf'][i], w=d['w irf'][i], t=p.t)
            x_m = np.array(self.params_dict2vec(tmp, separate_spacecraft=False))

        # Измерения с поправкой на угловой коэффициент усиления G (signal_rate)
        '''dr_cf, dr_ff, notes1 = measure_antennas_power(c=c, f=f, v=v, p=p, get_ready_measurements=True)
        z_ = np.array(dr_cf + dr_ff)
        # Поправка на G ниже
        signal_rate = np.ones(z_len)
        if v.NAVIGATION_ANGLES:
            signal_rate_cf, signal_rate_ff, notes2 = measure_antennas_power(c=c, f=f, v=v, p=p,
                                                                            get_signal_rates=True, j=j, x_m=x_m)
            signal_rate = np.array(signal_rate_cf + signal_rate_ff)
        z_ *= np.sqrt(signal_rate)'''
        z_ = v.MEASURES_VECTOR
        notes1 = v.MEASURES_VECTOR_NOTES
        # print(f"Измерения реальные:  {z_}")

        # Измерения согласно модели
        z_model, notes3 = measure_antennas_power(c=c, f=f, v=v, p=p, j=j, estimated_params=x_m)
        # print(f"Измерения модельные: {z_model}")
        # print(f"Положения модельные: {x_m}")

        # ZModel&RealDifference сейчас не относится к конкретным КА. Надо ли?
        # for i_c in range(c.n):
        #     for i_f in range(f.n):
        p.record.loc[p.iter, f'ZModel&RealDifference'] = np.abs(z_model - z_).mean()

        # Расчёт матриц
        Q_tilda = self.Phi @ self.D @ self.Q @ self.D.T @ self.Phi.T  # * self.v.dT  # nt_nt
        P_m = self.Phi @ self.P @ self.Phi.T + Q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        # >>>>>>>>>>>> Этап коррекции <<<<<<<<<<<<
        H, _, _, notesH = \
            h_matrix(c_ant=v.N_ANTENNA_C, f_ant=v.N_ANTENNA_F, fn=f.n, cn=c.n, r_f=f.r_orf, r_c=c.r_orf,
                     angles_navigation=v.NAVIGATION_ANGLES, multy_antenna_send=v.MULTI_ANTENNA_SEND,
                     multy_antenna_take=v.MULTI_ANTENNA_TAKE, w_0=v.W_ORB, t=p.t,
                     q_f=d['q-3 irf'], q_c=[c.q[i][1:4] for i in range(c.n)])

        R = np.eye(z_len) * v.KALMAN_COEF['r']  # n + n(n-1)/2
        my_print(f"P_m: {P_m.shape}, H: {R.shape}, P_m: {R.shape}", if_print=p.iter == 1)
        k_ = P_m @ H.T @ np.linalg.inv(H @ P_m @ H.T + R)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.P = (np.eye(j * f.n) - k_ @ H) @ P_m
        raw_estimation_params = np.array((np.matrix(x_m).T + k_ @ (z_ - z_model).T).T)[0]

        # >>>>>>>>>>>> Запись <<<<<<<<<<<<
        for i in range(f.n):
            tmp = raw_estimation_params[(0 + i) * j: (1 + i) * j]
            if v.SHAMANISM["KalmanQuaternionNormalize"] and v.NAVIGATION_ANGLES:
                tmp2 = vec2quat(tmp[3:6])
                tmp[3:6] = (tmp2 / np.linalg.norm(tmp2))[1:4]
            if v.SHAMANISM["KalmanSpinLimit"][0] and \
                    np.linalg.norm(tmp[9:12]) > v.SHAMANISM["KalmanSpinLimit"][1]:
                tmp[9:12] = tmp[9:12] / np.linalg.norm(tmp[9:12]) * v.SHAMANISM["KalmanSpinLimit"][1]
            if v.SHAMANISM["KalmanVelocityLimit"][0] and \
                    np.linalg.norm(tmp[6:9]) > v.SHAMANISM["KalmanVelocityLimit"][1]:
                tmp[6:9] = tmp[6:9] / np.linalg.norm(tmp[6:9]) * v.SHAMANISM["KalmanSpinLimit"][1]
            if v.SHAMANISM["KalmanPositionLimit"][0] and \
                    np.linalg.norm(tmp[0:3]) > v.SHAMANISM["KalmanPositionLimit"][1]:
                tmp[0:3] = tmp[0:3] / np.linalg.norm(tmp[0:3]) * v.SHAMANISM["KalmanPositionLimit"][1]
            self.estimation_params[i] = tmp

        # Запись и отображение
        if p.iter == 1:
            my_print(f"R-notes: {notes1}", color='y')
            my_print(f"Длина длин: {len(z_)}", color='r', if_print=v.IF_TEST_PRINT)
            my_print(f"M-notes: {notes3}", color='y')
            my_print(f"Длина модельных длин: {len(z_model)}", color='r', if_print=v.IF_TEST_PRINT)
            my_print(f"H-notes: {notesH}", color='y')
            with open("kiamfemto/data/measures_vector_notes_last.txt", "w") as f:
                f.write("# Рассчитано в PyCharm\n# Параметры: {rel} {N_1} {N_2} {i_1} {i_2} {send_len} {take_len}\n")
                f.write(f"# Параметр CUBESAT_AMOUNT {v.CUBESAT_AMOUNT}\n")
                f.write(f"# Параметр CHIPSAT_AMOUNT {v.CHIPSAT_AMOUNT}\n")
                f.write(f"# Параметр NAVIGATION_ANGLES {v.NAVIGATION_ANGLES}\n")
                f.write(f"# Параметр N_ANTENNA_C {v.N_ANTENNA_C}\n")
                f.write(f"# Параметр N_ANTENNA_F {v.N_ANTENNA_F}\n")
                f.write(f"# Параметр MULTI_ANTENNA_TAKE {v.MULTI_ANTENNA_TAKE}\n")
                f.write(f"# Параметр MULTI_ANTENNA_SEND {v.MULTI_ANTENNA_SEND}\n")
                for j in range(z_len):
                    f.write(f"{v.MEASURES_VECTOR_NOTES[j]}\n")
            my_print(f"P-: {P_m.shape}, H.T: {H.T.shape}, H: {H.shape}, R: {R.shape}", color='g')
            my_print(f"x-: {np.matrix(x_m).shape}, K: {k_.shape}, z: {(z_ - z_model).shape}", color='g')
            my_print(f"estimation_params: {raw_estimation_params.shape}", color='g')

def navigate(k: KalmanFilter):
    k.calc()  # Пока что при любом OPERATING_MODES (если все КА включены и не потеряны)

# >>>>>>>>>>>> Control <<<<<<<<<<<<
