"""Функции для моделирования динамики КА
ПЕРЕДЕЛАТЬ v_->v, v->u"""
from datetime import datetime

from primary_info import *
from gnc_systems import *
from my_plot import *

# >>>>>>>>>>>> Задание движения Хилла-Клохесси-Уилтшира <<<<<<<<<<<<
def get_c_hkw(r: Union[list, np.ndarray], v: Union[list, np.ndarray], w: float) -> list:
    """Возвращает константы C[0]..C[5] движения Хилла-Клохесси-Уилтштира"""
    return [2*r[2] + v[0]/w,
            v[2]/w,
            -3*r[2] - 2*v[0]/w,
            r[0] - 2*v[2]/w,
            v[1]/w,
            r[1]]

def r_hkw(C: Union[list, np.ndarray], w: float, t: float) -> np.ndarray:
    """Возвращает вектор координат в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    return np.array([-3*C[0]*w*t + 2*C[1]*np.cos(w*t) - 2*C[2]*np.sin(w*t) + C[3],
                     C[5]*np.cos(w*t) + C[4]*np.sin(w*t),
                     2*C[0] + C[2]*np.cos(w*t) + C[1]*np.sin(w*t)])

def v_hkw(C: Union[list, np.ndarray], w: float, t: float) -> np.ndarray:
    """Возвращает вектор скоростей в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    return np.array([-3 * C[0] * w - 2 * w * C[1] * np.sin(w * t) - 2 * w * C[2] * np.cos(w * t),
                     w * C[4] * np.cos(w * t) - w * C[5] * np.sin(w * t),
                     -w * C[2] * np.sin(w * t) + w * C[1] * np.cos(w * t)])

def get_rand_c(w: float, r_spread: float = 100, v_spread: float = 0.01, if_no: bool = False,
               if_quaternion: bool = False) -> list:
    """(quaternion or quaternion_but_i_dont_give_a_fuck)"""
    x, y, z = np.random.uniform(-r_spread, r_spread, 3)
    vx, vy, vz = np.random.uniform(-v_spread, v_spread, 3)
    a, b, g = np.random.uniform(-100, 100, 3)
    if if_quaternion and not if_no:
        return [2 * z + vx / w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y, a, b, g]
    return [2 * z + vx / w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y]

# >>>>>>>>>>>> Поступательное движение, интегрирование <<<<<<<<<<<<
def get_atm_params(v: Variables, h: float, atm_model: str = None) -> tuple:
    """NASA модель: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html (якорная точка: 25 км)
    ПНБО модель: https://www.energia.ru/ktt/archive/2022/04-2022/101-111.pdf (120-600 км)
    COESA62, COESA76 модели: библиотека poliastro
    :param v: Объект класса Variables
    :param h: Высота
    :param atm_model: Модель атмосферы (необязательный параметр)
    :return: (ρ, T, P): плотность, температура, давление (Внимание! Для ПНБО (ρ, None, None)"""
    atm_model = v.ATMOSPHERE_MODEL if atm_model is None else atm_model
    rho, T, p = None, None, None
    if atm_model == 'NASA':
        if h > 25e3:
            T = -131.21 + 0.00299 * h
            p = 2.488 * ((T + 273.1) / 216.6) ** -11.388
        elif h > 11e3:
            T = -56.46
            p = 22.65 * np.exp(1.73 - 0.000157 * h)
        else:
            T = 15.04 - 0.00649 * h
            p = 101.29 * (T + 273.1) / 288.08
        rho = p / (0.2869 * (T + 273.1))
    if atm_model == 'ПНБО':
        A = 2.519e-10
        H_m = 200e3
        H_0 = 290e3
        K_0 = 0.26
        a_1 = 100e3
        a_2 = 141.13e3
        n_0 = 6.34
        n_01 = 4.25
        n_02 = 4.37
        if h > 290e3:
            n = n_0 + K_0 * ((h - H_0) / a_1) ** n_01 - ((h - H_0) / a_2) ** n_02
        else:
            n = n_0 + K_0 * ((H_0 - h) / a_1) ** n_01
        rho = A * (H_m / h) ** n
    if atm_model in ['COESA62', 'COESA76']:
        from astropy import units as u
        from poliastro.earth.atmosphere import COESA62, COESA76
        coesa = COESA62() if atm_model == 'COESA62' else COESA76()
        T, p, rho = coesa.properties(h * u.m)
        T = T.value
        p = p.value
        rho = rho.value
    return rho, T, p

def get_geopotential_acceleration(v_: Variables, r: Union[list, np.ndarray], v: Union[list, np.ndarray]) -> np.ndarray:
    """Возвращает ускорение КА от притяжения Земли.
    Внимание! При глобальном параметре DYNAMIC_MODEL='Clohessy-Wiltshire' возвращает ускорение в ОСК.
    Иначе возвращает ускорение в ИСК.
    :param v: Скорость КА. Внимание! Не путать с Variables!
    :param r: Положение КА
    :param v_: Переменная класса Variables. Внимание! Не путать со скоростью!"""
    if 'hkw' in v_.SOLVER:
        return np.array([-2 * v_.W_ORB * v[2],
                         -v_.W_ORB ** 2 * r[1],
                         2 * v_.W_ORB * v[0] + 3 * v_.W_ORB ** 2 * r[2]])
    return v_.MU * np.array(r) / np.linalg.norm(r)**3

def get_aero_drag_acceleration(v_: Variables, obj: Union[CubeSat, FemtoSat], i: int,
                               r: Union[list, np.ndarray], v: Union[list, np.ndarray]):
    """Возвращает ускорение КА от сопротивления атмосферы.
    Внимание! При глобальном параметре DYNAMIC_MODEL='Clohessy-Wiltshire' возвращает ускорение в ОСК.
    Иначе возвращает ускорение в ИСК."""
    S = quart2dcm(obj.q[i])  # Для новых кватернионов - неверно!
    cos_alpha = clip((np.trace(S) - 1) / 2, -1, 1)
    # alpha = 180 / np.pi * np.arccos(cos_alpha)
    rho = get_atm_params(v=v_, h=obj.r_orf[i][2] + v_.ORBIT_RADIUS - v_.EARTH_RADIUS)[0]
    if obj.name == "CubeSat":
        c_resist = 1.05
        square = obj.size[0] * obj.size[1]
    else:
        c_resist = 1.17
        square = obj.size[0] * obj.size[1] * abs(cos_alpha)

    if 'hkw' in v_.SOLVER:
        v_real = v + np.array([v_.V_ORB, 0, 0])
        rho = get_atm_params(v=v_, h=r[2] + v_.ORBIT_RADIUS - v_.EARTH_RADIUS)[0]
        return - v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / obj.mass

def get_full_acceleration(v_: Variables, obj: Union[CubeSat, FemtoSat], i: int,
                          r: Union[float, np.ndarray], v: Union[float, np.ndarray]) -> np.ndarray:
    """Возвращает вектор силы в ОСК, принимает параметры в ОСК\n
    square - площадь S для аэродинамики
    c_resist - Cf для аэродинаммики"""
    if 'hkw' in v_.SOLVER:
        force = get_geopotential_acceleration(v_=v_, r=r, v=v)
        if v_.DYNAMIC_MODEL['aero drag']:
            force += get_aero_drag_acceleration(v_=v_, r=r, v=v, obj=obj, i=i)
        return force

def rk4_translate(v_: Variables, obj: Union[CubeSat, FemtoSat], i: int, dt: float = None, r=None, v=None) -> tuple:
    dt = v_.dT if dt is None else dt

    def rv_right_part(rv1, a1):
        return np.array([rv1[3], rv1[4], rv1[5], a1[0], a1[1], a1[2]])
    r = obj.r_orf[i] if r is None else r
    v = obj.v_orf[i] if v is None else v
    a = get_full_acceleration(v_=v_, obj=obj, i=i, r=r, v=v)

    rv = np.append(r, v)
    k1 = rv_right_part(rv, a)
    k2 = rv_right_part(rv + k1 * dt / 2, a)
    k3 = rv_right_part(rv + k2 * dt / 2, a)
    k4 = rv_right_part(rv + k3 * dt, a)
    rv = dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return rv[0:3] + r, rv[3:6] + v


# >>>>>>>>>>>> Вращательное движение, интегрирование <<<<<<<<<<<<
def get_torque(v_: Variables, obj: Union[CubeSat, FemtoSat], q: Union[list, np.ndarray], w: np.ndarray) -> np.ndarray:
    """Возвращает вектор углового УСКОРЕНИЯ. Удачи!"""
    return np.zeros(3)

def rk4_attitude(v_: Variables, obj: Union[CubeSat, FemtoSat], t: float, i: int, dt: float = None, q=None, w=None):
    """Господи, где здесь производная, где дифференциал? Пожалуйста, дрогой я, дай мне знак! ДОрОжНЫй!!!"""
    dt = v_.dT if dt is None else dt

    U, S, A, R_orb = get_matrices(v=v_, t=t, obj=obj, n=i)

    def lw_right_part(qw_, e_):
        q_, w_ = qw_[0:4], qw_[4:7]
        dq = 1 / 2 * q_dot(q_, [0, w_[0], w_[1], w_[2]])
        # J = A.T @ obj.J @ A
        J = obj.J
        dw = - (np.linalg.inv(J) @ (my_cross(w_, J @ w_)))  # + e_
        return np.append(dq, dw)
    q = obj.q[i] if q is None else q
    w = obj.w_irf[i] if w is None else w
    e = get_torque(v_=v_, obj=obj, q=q, w=w)

    a = len(q)
    q4 = q if a == 4 else vec2quat(q)
    qw = np.append(q4, w)
    k1 = lw_right_part(qw, e)
    k2 = lw_right_part(qw + k1 * dt / 2, e)
    k3 = lw_right_part(qw + k2 * dt / 2, e)
    k4 = lw_right_part(qw + k3 * dt, e)
    qw = dt / 6 * (k1 + 2*k2 + 2*k3 + k4)
    # q_anw = q_dot(qw[0:4]/np.linalg.norm(qw[0:4]), q4)
    q_anw = q4 + qw[0:4]
    q_anw /= np.linalg.norm(q_anw)
    return q_anw[4 - a:4], qw[a:a + 3] + w


# >>>>>>>>>>>> Перевод между системами координат <<<<<<<<<<<<
def get_matrices(v: Variables, t: float, obj: Apparatus, n: int, first_init: bool = False):
    """Функция возвращает матрицы поворота.
    Инициализируется в dymancis.py, используется в spacecrafts, dynamics"""
    E = t * v.W_ORB  # Эксцентрическая аномалия
    f = 2 * np.arctan(np.sqrt((1 + v.ECCENTRICITY) / (1 - v.ECCENTRICITY)) * np.tan(E / 2))  # Истинная аномалия
    A = quart2dcm(obj.q[n])
    if 'hkw' in v.SOLVER or first_init:
        U = np.array([[0, 1, 0],  # Поворот к экваториальной плоскости
                      [0, 0, 1],
                      [1, 0, 0]]) @ \
            np.array([[np.cos(f), np.sin(f), 0],  # Разница между истинной аномалией и местной
                      [-np.sin(f), np.cos(f), 0],
                      [0, 0, 1]]) @ \
            np.array([[1, 0, 0],
                      [0, np.cos(deg2rad(v.INCLINATION)), np.sin(deg2rad(v.INCLINATION))],  # Поворот к плоскости орбиты
                      [0, -np.sin(deg2rad(v.INCLINATION)), np.cos(deg2rad(v.INCLINATION))]])
        translation = v.P / (1 + v.ECCENTRICITY * np.cos(f))
    else:
        e_z = vec2unit(v.ANCHOR.r_irf[0])
        e_x = vec2unit(v.ANCHOR.v_irf[0])
        e_y = vec2unit(np.cross(e_z, e_x))
        e_x = vec2unit(np.cross(e_y, e_z))
        U = np.array([e_x, e_y, e_z])
        translation = np.linalg.norm(v.ANCHOR.r_irf[0])
    S = A @ U.T
    # R_orb = U.T @ np.array([0, 0, v.ORBIT_RADIUS])
    R_orb = U.T @ np.array([0, 0, translation])
    return U, S, A, R_orb

def i_o(a, v: Variables, U: np.ndarray, vec_type: str):
    """Инерциальная -> Орбитальная"""
    a_np = np.array(a)
    if len(a_np.shape) == 1 and vec_type == "r":
        return U @ a_np - np.array([0, 0, v.ORBIT_RADIUS])
    if len(a_np.shape) == 1 and vec_type == "v":
        return U @ a_np - np.array([v.V_ORB, 0, 0])
    if len(a_np.shape) == 1 and vec_type == "w":
        return U @ (a_np - v.W_ORB_VEC_IRF)
    if len(a_np.shape) == 2:
        return U @ a_np @ U.T
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v, w]")

def o_i(a, v: Variables, U: np.ndarray, vec_type: str):
    """Орбитальная -> Инерциальная"""
    a_np = np.array(a)
    if len(a_np.shape) == 1 and vec_type == "r":
        return U.T @ (a_np + np.array([0, 0, v.ORBIT_RADIUS]))
    if len(a_np.shape) == 1 and vec_type == "v":
        return U.T @ (a_np + np.array([v.V_ORB, 0, 0]))
    if len(a_np.shape) == 1 and vec_type == "w":
        return U.T @ a_np + v.W_ORB_VEC_IRF
    if len(a_np.shape) == 2:
        return U.T @ a_np @ U
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v, w]")


# >>>>>>>>>>>> Класс динамики кубсатов и чипсатов <<<<<<<<<<<<
class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat, a: Anchor, v: Variables):
        from pandas import DataFrame

        # Неизменные параметры
        self.t = 0.
        self.iter = 0
        self.v = v
        self.c = c
        self.f = f
        self.a = a
        self.spacecrafts_cd = [self.c, self.f]
        self.spacecrafts_all = [self.a, self.c, self.f]
        self.time_begin = datetime.now()

        # Инициализация фильтра
        self.k = KalmanFilter(f=f, c=c, p=self)

        # Инициализация траектории kiam-astro
        self.jd0, self.tr = None, None
        if 'kiamastro' in self.v.SOLVER:
            self.kiam_init()

        # Запись параметров
        self.record = DataFrame()
        self.do_report()  # Продублировать в конце time_step()
        self.record = self.record.astype({'i': 'int32', 'FemtoSat n': 'int32', 'CubeSat n': 'int32'})

    def kiam_init(self):
        from kiam_astro import kiam
        from kiam_astro.trajectory import Trajectory
        self.jd0 = kiam.juliandate(2024, 1, 1, 0, 0, 0)  # (год, месяц, день, чч, мм, сс)
        self.tr = [[Trajectory(initial_state=np.zeros(6), initial_time=0, initial_jd=self.jd0, variables='rv',
                               system='gcrs', units_name='earth') for _ in range(obj.n)]
                   for obj in [self.c, self.f, self.a]]
        for j, obj in enumerate(self.spacecrafts_all):
            for i in range(obj.n):
                s0 = np.append(obj.r_irf[i] / (kiam.units('earth')['DistUnit'] * 1e3),
                               obj.v_irf[i] / (kiam.units('earth')['VelUnit'] * 1e3))
                self.tr[j][i] = Trajectory(initial_state=s0, initial_time=0, initial_jd=self.jd0, variables='rv',
                                           system='gcrs', units_name='earth')
                self.tr[i][j].set_model(variables='rv', model_type='nbp', primary='earth',
                                        sources_list=[] + ['j2'] if self.v.DYNAMIC_MODEL['j2'] else [] +
                                                     ['atm'] if self.v.DYNAMIC_MODEL['aero drag'] else [])
                self.tr[i][j].model['data']['jd_zero'] = self.jd0
                self.tr[i][j].model['data']['mass'] = self.f.mass
                self.tr[i][j].model['data']['area'] = self.f.size[0] * self.f.size[1]
                self.tr[i][j].model['data']['order'] = 0  # order of the Moon's gravity field
                self.tr[i][j].propagate(tof=self.v.TIME/self.v.SEC_IN_RAD, npoints=int(self.v.TIME//self.v.dT))
        my_print(f"kiam-astro Time: {self.tr[0][0].times[-1]} ({self.tr[0][0].times[-1]/2/np.pi} оборотов)\n"
                 f"kiam-astro Points: {self.v.TIME/self.v.dT}", if_print=self.v.IF_ANY_PRINT)
        if self.v.IF_ANY_SHOW:
            self.tr[0][0].show(variables='3d', language='rus')
            self.tr[1][0].show(variables='3d', language='rus')

    # Шаг по времени
    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.v.dT

        # Движение системы
        for j, obj in enumerate(self.spacecrafts_all):
            for i in range(obj.n):
                # Вращательное движение
                if obj != self.a:
                    obj.q[i], obj.w_irf[i] = rk4_attitude(v_=self.v, obj=obj, i=i, t=self.t)

                # Поступательное движение
                if 'rk4' in self.v.SOLVER:
                    if np.any(self.v.DYNAMIC_MODEL.values()):  # Если J2 или aero drag
                        obj.r_orf[i], obj.v_orf[i] = rk4_translate(v_=self.v, obj=obj, i=i)
                    else:
                        obj.r_orf[i] = r_hkw(obj.c_hkw[i], self.v.W_ORB, self.t)
                        obj.v_orf[i] = v_hkw(obj.c_hkw[i], self.v.W_ORB, self.t)

                    U, _, _, _ = get_matrices(v=self.v, t=self.t, obj=obj, n=i)
                    obj.r_irf[i] = o_i(v=self.v, a=obj.r_orf[i], U=U, vec_type='r')
                    obj.v_irf[i] = o_i(v=self.v, a=obj.v_orf[i], U=U, vec_type='v')
                elif 'kiamastro' in self.v.SOLVER:
                    from kiam_astro import kiam
                    obj.r_irf[i] = np.array([self.tr[j][i].states[ii][self.iter - 1]
                                             for ii in range(3)]) * kiam.units('earth')['DistUnit'] * 1e3
                    obj.v_irf[i] = np.array([self.tr[j][i].states[ii + 3][self.iter - 1]
                                             for ii in range(3)]) * kiam.units('earth')['VelUnit'] * 1e3
                    tr_time = self.tr[j][i].times[self.iter-1] * self.v.SEC_IN_RAD
                    U, _, _, _ = get_matrices(v=self.v, t=tr_time, obj=obj, n=i)
                    obj.r_orf[i] = i_o(v=self.v, a=obj.r_irf[i], U=U, vec_type='r')
                    obj.v_orf[i] = i_o(v=self.v, a=obj.v_irf[i], U=U, vec_type='v')
                else:
                    raise ValueError(f"Поменяй солвер! SOLVER={self.v.SOLVER}, а должен быть среди {self.v.SOLVERS}!")

                obj.w_orf[i] = i_o(v=self.v, a=obj.w_irf[i], U=U, vec_type='w')

        # Комплекс первичной информации
        measure_antennas_power(c=self.c, f=self.f, v=self.v, noise=np.sqrt(self.v.KALMAN_COEF['r']), produce=True,
                               p=self)
        measure_magnetic_field(c=self.c, f=self.f, v=self.v, noise=np.sqrt(self.v.KALMAN_COEF['r']))

        # Изменение режимов работы
        guidance(v=self.v, c=self.c, f=self.f, earth_turn=self.t * self.v.W_ORB / 2 / np.pi)

        # Запись параметров
        self.do_report()

        # Навигация чипсатов
        if self.v.IF_NAVIGATION:
            navigate(k=self.k)

    def do_report(self):
        i_t = self.iter
        d = self.record
        d.loc[i_t, f'i'] = self.iter
        d.loc[i_t, f't'] = self.t
        for obj in self.spacecrafts_cd:
            d.loc[i_t, f'{obj.name} n'] = obj.n
            for i_n in range(obj.n):
                for v in ['r', 'q', 'v', 'w']:
                    tmp = {'r': [obj.r_irf[i_n], obj.r_orf[i_n]],
                           'v': [obj.v_irf[i_n], obj.v_orf[i_n]],
                           'q': [obj.q[i_n][1:4]],
                           'w': [obj.w_irf[i_n], obj.w_orf[i_n]]}[v]
                    for i_fr, frame in enumerate(['irf', 'orf'] if v != 'q' else ['irf']):  # Кватернионы только в ИСК
                        for i_r, c in enumerate('xyz'):
                            d.loc[i_t, f'{obj.name} {v} {c} {frame} {i_n}'] = tmp[i_fr][i_r]

        for obj in [self.f]:
            for i_n in range(obj.n):
                if obj.operating_mode[i_n] != "lost":  # Иначе заполняется Null (в plot в self.v.NO_LINE_FLAG)
                    r_orf_estimation = self.k.get_estimation(i_f=i_n, v='r orf')
                    w_irf_estimation = self.k.get_estimation(i_f=i_n, v='w irf')
                    q_irf_estimation = self.k.get_estimation(i_f=i_n, v='q-3 irf')
                    r_orf = self.f.r_orf[i_n]
                    w_irf = self.f.w_irf[i_n]
                    q_irf = self.f.q[i_n][1:4]

                    if self.v.NAVIGATION_ANGLES:
                        U, _, _, _ = get_matrices(v=self.v, t=self.t, obj=obj, n=i_n)
                        w_orf = i_o(a=w_irf, v=self.v, vec_type='w', U=U)
                        w_orf_estimation = i_o(a=w_irf_estimation, v=self.v, vec_type='w', U=U)

                    d.loc[i_t, f'{obj.name} KalmanPosEstimation r {i_n}'] = np.linalg.norm(r_orf_estimation)
                    d.loc[i_t, f'{obj.name} KalmanPosError r {i_n}'] = np.linalg.norm(r_orf_estimation - r_orf)
                    if self.v.NAVIGATION_ANGLES:
                        d.loc[i_t, f'{obj.name} KalmanSpinError w {i_n}'] = np.linalg.norm(w_irf_estimation - w_irf)
                        d.loc[i_t, f'{obj.name} KalmanQuatError q {i_n}'] = np.linalg.norm(q_irf_estimation - q_irf)
                    for i_r, c in enumerate('xyz'):
                        d.loc[i_t, f'{obj.name} KalmanPosEstimation {c} {i_n}'] = r_orf_estimation[i_r]
                        d.loc[i_t, f'{obj.name} KalmanPosError {c} {i_n}'] = r_orf_estimation[i_r] - r_orf[i_r]
                        if self.v.NAVIGATION_ANGLES:
                            d.loc[i_t, f'{obj.name} RealSpin IRF {c} {i_n}'] = w_irf[i_r]
                            d.loc[i_t, f'{obj.name} RealSpin ORF {c} {i_n}'] = w_orf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinEstimation IRF {c} {i_n}'] = w_irf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinEstimation ORF {c} {i_n}'] = w_orf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinError IRF {c} {i_n}'] = w_irf_estimation[i_r] - w_irf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanSpinError ORF {c} {i_n}'] = w_orf_estimation[i_r] - w_orf[i_r]
                            d.loc[i_t, f'{obj.name} RealQuat {c} {i_n}'] = q_irf[i_r]
                            d.loc[i_t, f'{obj.name} KalmanQuatEstimation {c} {i_n}'] = q_irf_estimation[i_r]
                            d.loc[i_t, f'{obj.name} KalmanQuatError {c} {i_n}'] = q_irf_estimation[i_r] - q_irf[i_r]
