"""Ёбаный пиздец блять в каком же я ахуе просто ебааать"""
from srs.kiamfemtosat.spacecrafts import *
from srs.kiamfemtosat.data.observability_mapping_partial_derivatives import *

# >>>>>>>>>>>> Guidance <<<<<<<<<<<<
def guidance(c: CubeSat, f: FemtoSat, v: Variables, earth_turn: float) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            if obj.operating_modes[i] == v.OPERATING_MODES_CHANGE[1]:  # Отсутствие аккумулятора на чипсате
                if earth_turn % 1 < 0.5 and obj.operating_mode[i] == v.OPERATING_MODES[-1]:
                    obj.operating_mode[i] = v.OPERATING_MODES[0]
                if earth_turn % 1 > 0.5 and obj.operating_mode[i] != v.OPERATING_MODES[-1]:
                    obj.operating_mode[i] = v.OPERATING_MODES[-1]

# >>>>>>>>>>>> Navigation <<<<<<<<<<<<
class KalmanFilter:
    def __init__(self, f: FemtoSat, c: CubeSat, p: any):
        # Общие параметры
        self.f = f  # Фемтоспутники
        self.c = c  # Кубсаты
        self.p = p  # Динамическая модель
        self.v = p.v
        self.j = len(f.rv_orf_calc[0]) if self.v.NAVIGATION_ANGLES else 6  # Вектор состояния 1 чипсата
        self.sequence_under_diagonal = flatten([[[i, j] for i in range(j)] for j in range(self.f.n)])
        self.sigmas, self.real_sigmas = [[[] for _ in range(self.j * self.f.n)] for _ in range(2)]  # rename -> errors

        # Матрицы фильтра в начальный момент времени
        if not self.v.NAVIGATION_ANGLES:
            self.r_orf_estimation = np.array([np.append(f.rv_orf_calc[i][0:3],
                                                        f.rv_orf_calc[i][7:10]) for i in range(f.n)])
            self.Phi = np.array([[0, 0, 0, 1, 0, 0],
                                 [0, 0, 0, 0, 1, 0],
                                 [0, 0, 0, 0, 0, 1],
                                 [0, 0, 0, 0, 0, -2*self.v.W_ORB],
                                 [0, - self.v.W_ORB**2, 0, 0, 0, 0],
                                 [0, 0, 3*self.v.W_ORB**2,  2*self.v.W_ORB, 0, 0]]) * self.v.dT + np.eye(6)  # !!!!!!!!!
            self.D = np.vstack([np.zeros((3, 3)), np.eye(3)])
            self.P = [np.diag([self.v.KALMAN_COEF['p'][0]] * 3 + [self.v.KALMAN_COEF['p'][1]] * 3) for _ in range(f.n)]
            self.Q = np.diag([self.v.KALMAN_COEF['q'][0]] * 3)
        else:  # Вектор состояния содержит угловые переменные и скорости
            self.r_orf_estimation = np.array(f.rv_orf_calc)
            self.Phi = np.array([[0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, -2*self.v.W_ORB, 0, 0, 0],
                                 [0, -self.v.W_ORB**2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 3*self.v.W_ORB**2, 0, 0, 0, 0, 2*self.v.W_ORB, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]) * self.v.dT + np.eye(self.j)  # !!!!!!!!!!!!!
            self.P = [np.diag([self.v.KALMAN_COEF['p'][0]] * 3 + [self.v.KALMAN_COEF['p'][2]] * 4 +
                              [self.v.KALMAN_COEF['p'][1]] * 3 + [self.v.KALMAN_COEF['p'][3]] * 3) for _ in range(f.n)]
            self.D = np.vstack([np.zeros((7, 6)), np.eye(6)])
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
                              [np.zeros([self.j, tmp])] * (self.f.n - i - 1) for i in range(self.f.n)])

        # Вывод
        if self.v.IF_ANY_PRINT:
            my_print(f"Матрицы Ф:{self.Phi.shape}, Q:{self.Q.shape}, P:{self.P.shape}, D:{self.D.shape}", color='g')

    def calc(self) -> None:  # Считается, что NAVIGATION_BY_ALL = True
        from srs.kiamfemtosat.primary_info import measure_antennas_power
        from srs.kiamfemtosat.dynamics import rk4_translate, rk4_attitude

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
        # Моделирование орбитального движения на dT -> вектор состояния x_m
        if self.v.DYNAMIC_MODEL['aero drag']:  # ПЕРЕДЕЛАТЬ ВЕСЬ РАСЧЁТ (вроде всё хорошо, примечание можно убрать)
            rv_m = [rk4_translate(v_=self.v, obj=self.f, i=i, r=self.r_orf_estimation[i][0:3],
                                  v=self.r_orf_estimation[i][7:10] if self.v.NAVIGATION_ANGLES else
                                  self.r_orf_estimation[i][3:6]) for i in range(self.f.n)]
            tmp = int(self.v.NAVIGATION_ANGLES)
            x_m = np.array(flatten([np.append(np.append(np.append(rv_m[i][0], [0]*4*tmp), rv_m[i][1]), [0]*3*tmp)
                                    for i in range(self.f.n)]))
        else:
            x_m = np.array(self.Phi @ np.array(flatten(self.r_orf_estimation)))[0]  # flatten -> [1..j] + ... (f.n раз)

        # Моделирование углового движения на dT -> вектор состояния x_m
        if self.v.NAVIGATION_ANGLES:
            for i in range(self.f.n):
                qw_m = rk4_attitude(v_=self.v, obj=self.f, i=i,
                                    q=self.r_orf_estimation[i][3:7], w=self.r_orf_estimation[i][10:13])
                x_m[i*self.j + 3:i*self.j + 7], x_m[i*self.j + 10:i*self.j + 13] = (qw_m[0], qw_m[1])

        # Измерения с поправкой на угловой коэффициент усиления G (signal_rate)
        dr_cf, dr_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v, get_ready_measurements=True)
        my_print(f"Длина длин 1: {len(dr_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        my_print(f"Длина длин 2: {len(dr_cf + dr_ff)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        z_ = np.array(dr_cf + dr_ff)
        ''''for i_c in range(self.c.n):
            for i_f in range(self.f.n):
                cf_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_take=True)) * \
                         len(get_gain(v=self.v, obj=self.c, r=np.ones(3), if_send=True))
                z_ = np.append(z_, self.v.MEASURES_VECTOR[count:count+cf_len])
                count += cf_len
                fc_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_send=True)) * \
                         len(get_gain(v=self.v, obj=self.c, r=np.ones(3), if_take=True))
                z_ = np.append(z_, self.v.MEASURES_VECTOR[count:count+fc_len])
                count += fc_len
        for i_f1 in range(self.f.n):
            for i_f2 in range(i_f1):
                if i_f1 != i_f2:
                    ff_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_take=True)) * \
                             len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_send=True))
                    z_ = np.append(z_, self.v.MEASURES_VECTOR[count:count+ff_len])
                    count += ff_len
                    ff_len = len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_send=True)) * \
                             len(get_gain(v=self.v, obj=self.f, r=np.ones(3), if_take=True))
                    z_ = np.append(z_, self.v.MEASURES_VECTOR[count:count+ff_len])
                    count += ff_len'''

        if self.v.NAVIGATION_ANGLES:
            signal_rate_cf, signal_rate_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v,
                                                                           get_signal_rates=True, j=self.j, x_m=x_m)
            my_print(f"Длина G 1: {len(signal_rate_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
            my_print(f"Длина G 2: {len(signal_rate_cf + signal_rate_ff)}", color='r',
                     if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
            signal_rate = np.array(signal_rate_cf + signal_rate_ff)
        else:
            signal_rate = np.array([1] * z_len)
            ''''for i_c in range(self.c.n):
                for i_f in range(self.f.n):
                    # Измерения чипсата сигналов от кубсата
                    tmp1 = get_gain(v=self.v, obj=self.f, if_take=True,
                                    r=quart2dcm(x_m[i_f*self.j+3:i_f*self.j+7]) @
                                      (x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c]))
                    tmp2 = get_gain(v=self.v, obj=self.c, if_send=True,
                                    r=quart2dcm(self.c.q[i_c]) @ (x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c]))
                    signal_rate = np.append(signal_rate, flatten([[tmp1[ii] * tmp2[jj] for ii in range(len(tmp1))]
                                                                                       for jj in range(len(tmp2))]))
                    # Измерения кубсата сигналов от чипсата
                    tmp1 = get_gain(v=self.v, obj=self.f, if_send=True,
                                    r=quart2dcm(x_m[i_f*self.j+3:i_f*self.j+7]) @
                                      (x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c]))
                    tmp2 = get_gain(v=self.v, obj=self.c, if_take=True,
                                    r=quart2dcm(self.c.q[i_c]) @ (x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c]))
                    signal_rate = np.append(signal_rate, flatten([[tmp1[ii] * tmp2[jj] for ii in range(len(tmp1))]
                                                                                       for jj in range(len(tmp2))]))
            for i_f1 in range(self.f.n):
                for i_f2 in range(i_f1):
                    if i_f1 != i_f2:
                        # Измерения чипсата 1 сигналов от чипсата 2
                        tmp1 = get_gain(v=self.v, obj=self.f, if_take=True,
                                        r=quart2dcm(x_m[i_f1*self.j+3:i_f1*self.j+7]) @
                                          (x_m[i_f1*self.j+0:i_f1*self.j+3] - x_m[i_f2*self.j+0:i_f2*self.j+3]))
                        tmp2 = get_gain(v=self.v, obj=self.f, if_send=True,
                                        r=quart2dcm(x_m[i_f2*self.j+3:i_f2*self.j+7]) @
                                          (x_m[i_f1*self.j+0:i_f1*self.j+3] - x_m[i_f2*self.j+0:i_f2*self.j+3]))
                        signal_rate = np.append(signal_rate, flatten([[tmp1[ii] * tmp2[jj] for ii in range(len(tmp1))]
                                                                                           for jj in range(len(tmp2))]))
                        # Измерения чипсата 2 сигналов от чипсата 1
                        tmp1 = get_gain(v=self.v, obj=self.f, if_send=True,
                                        r=quart2dcm(x_m[i_f1*self.j+3:i_f1*self.j+7]) @
                                          (x_m[i_f1*self.j+0:i_f1*self.j+3] - x_m[i_f2*self.j+0:i_f2*self.j+3]))
                        tmp2 = get_gain(v=self.v, obj=self.f, if_take=True,
                                        r=quart2dcm(x_m[i_f2*self.j+3:i_f2*self.j+7]) @
                                          (x_m[i_f1*self.j+0:i_f1*self.j+3] - x_m[i_f2*self.j+0:i_f2*self.j+3]))
                        signal_rate = np.append(signal_rate, flatten([[tmp2[ii] * tmp1[jj] for ii in range(len(tmp2))]
                                                                                           for jj in range(len(tmp1))]))
                                                                                           '''
        z_ *= np.sqrt(signal_rate)

        # Измерения согласно модели
        z_cf, z_ff, notes = measure_antennas_power(c=self.c, f=self.f, v=self.v, get_model_state=True, j=self.j,
                                                   x_m=x_m)
        my_print(f"Длина модельных длин 1: {len(z_cf)}", color='r', if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        my_print(f"Длина модельных длин 2: {len(z_cf + z_ff)}", color='r',
                 if_print=self.p.iter == 1 and self.v.IF_TEST_PRINT)
        z_model = np.array([z_cf + z_ff])
        '''
        for i_c in range(self.c.n):
            for i_f in range(self.f.n):
                # К чипсату от кубсата
                cf_len = len(get_gain(v=self.v, obj=self.f, r=np.random.random(3), if_take=True)) * \
                         len(get_gain(v=self.v, obj=self.c, r=np.random.random(3), if_send=True))
                z_model = np.append(z_model, [np.linalg.norm(x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c])
                                              for _ in range(cf_len)])
                # К кубсату от чипсата
                fc_len = len(get_gain(v=self.v, obj=self.f, r=np.random.random(3), if_send=True)) * \
                         len(get_gain(v=self.v, obj=self.c, r=np.random.random(3), if_take=True))
                z_model = np.append(z_model, [np.linalg.norm(x_m[i_f*self.j+0:i_f*self.j+3] - self.c.r_orf[i_c])
                                              for _ in range(fc_len)])
        for i_f1 in range(self.f.n):
            for i_f2 in range(i_f1):
                ff_len = len(get_gain(v=self.v, obj=self.f, r=np.random.random(3), if_take=True)) * \
                         len(get_gain(v=self.v, obj=self.f, r=np.random.random(3), if_send=True))
                z_model = np.append(z_model, [np.linalg.norm(x_m[i_f2*self.j+0:i_f2*self.j+3] -
                                                             x_m[i_f1*self.j+0:i_f1*self.j+3]) for _ in range(ff_len)])
                z_model = np.append(z_model, [np.linalg.norm(x_m[i_f2*self.j+0:i_f2*self.j+3] -
                                                             x_m[i_f1*self.j+0:i_f1*self.j+3]) for _ in range(ff_len)])
        '''        
        for i in range(self.f.n):  # ПРОВЕРИТЬ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tmp1 = get_gain(v=self.v, obj=self.f, r=np.random.random(3))
            self.f.z_difference[i] += [np.linalg.norm(z_model[i*2*len(tmp1): i*2*len(tmp1) + 1] -
                                                      z_[i*2*len(tmp1): i*2*len(tmp1) + 1])]

        # Расчёт матриц
        Q_tilda = self.Phi @ self.D @ self.Q @ self.D.T @ self.Phi.T * self.v.dT  # nt_nt
        P_m = self.Phi @ self.P @ self.Phi.T + Q_tilda  # nt_nt @ nt_nt @ nt_nt -> nt_nt

        H, _, _ = h_matrix(c_ant=self.v.N_ANTENNA_C, f_ant=self.v.N_ANTENNA_F, fn=self.f.n, cn=self.c.n,
                           angles_navigation=self.v.NAVIGATION_ANGLES, r_f=self.f.r_orf, r_c=self.c.r_orf,
                           multy_antenna_send=self.v.MULTI_ANTENNA_SEND, multy_antenna_take=self.v.MULTI_ANTENNA_TAKE,
                           w_0=self.v.W_ORB, t=self.p.t, q_f=self.f.q, q_c=self.c.q)
        if self.p.iter == 1:
            for j in range(z_len):
                print(self.v.MEASURES_VECTOR_NOTES[j])

        R = np.eye(z_len) * self.v.KALMAN_COEF['r']  # n + n(n-1)/2
        my_print(f"P-: {P_m.shape}, H.T: {H.T.shape}, H: {H.shape}, R: {R.shape}",
                 color='g', if_print=self.p.iter == 1)
        k_ = P_m @ H.T @ np.linalg.inv(H @ P_m @ H.T + R)  # nt_nt @ nt_lt @ l_l -> nt_l
        self.P = (np.eye(self.j * self.f.n) - k_ @ H) @ P_m
        my_print(f"x-: {np.matrix(x_m).shape}, K: {k_.shape}, z: {(z_ - z_model).shape}",
                 color='g', if_print=self.p.iter == 1)
        r_orf_estimation = np.matrix(x_m) + k_ @ (z_ - z_model).T
        for i in range(self.f.n):
            tmp = np.array(r_orf_estimation)[0][(0 + i)*self.j:(1 + i)*self.j]
            if self.v.SHAMANISM["KalmanQuaternionNormalize"] and self.v.NAVIGATION_ANGLES:
                tmp[3:7] = tmp[3:7] / np.linalg.norm(tmp[3:7])
            if self.v.SHAMANISM["KalmanSpinLimit"][0] and \
                    np.linalg.norm(tmp[10:13]) > self.v.SHAMANISM["KalmanSpinLimit"][1]:
                tmp[10:13] = tmp[10:13] / np.linalg.norm(tmp[10:13]) * self.v.SHAMANISM["KalmanSpinLimit"][1]
            self.r_orf_estimation[i] = tmp

def navigate(k: KalmanFilter):
    if k.v.NAVIGATION_BY_ALL:  # self.k.single_femto_filter:
        k.calc()  # Пока что при любом OPERATING_MODES (если весь рой выпал)
    else:
        for i in range(k.f.n):
            if k.f.operating_mode[i] != k.v.OPERATING_MODES[-1]:
                k.calc(i)
            else:
                k.f.z_difference[i] += [k.v.NO_LINE_FLAG]
            if k.c.gain_mode != k.v.GAIN_MODES[4]:
                for j_n in range(k.f.n):
                    for j_t in range(9 if k.v.NAVIGATION_ANGLES else 3):  # range(self.k.t)
                        tmp = np.append(np.append(k.f.r_orf[j_n], k.f.v_orf[j_n]), list(k.f.q[j_n][1:4]))
                        tmp = np.zeros(9) if k.f.operating_mode[j_n] != k.v.OPERATING_MODES[-1] else tmp
                        k.sigmas[j_n * k.j + j_t] += [np.sqrt(k.P[j_n][j_t][j_t]) * tmp[j_t]]

# >>>>>>>>>>>> Control <<<<<<<<<<<<
