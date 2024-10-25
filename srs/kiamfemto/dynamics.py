"""Функции для моделирования динамики КА
ПЕРЕДЕЛАТЬ v_->v, v->u"""
from datetime import datetime

from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

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

def rk4_attitude(v_: Variables, obj: Union[CubeSat, FemtoSat], i: int, dt: float = None, q=None, w=None):
    """Господи, где здесь производная, где дифференциал? Пожалуйста, дрогой я, дай мне знак! ДОрОжНЫй!!!"""
    dt = v_.dT if dt is None else dt

    def lw_right_part(qw_, e_):
        q_, w_ = qw_[0:4], qw_[4:7]
        dq = 1 / 2 * q_dot([0, w_[0], w_[1], w_[2]], q_)
        return np.append(dq, e_)
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
    qw = dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    tmp = (qw[0:4] + q4) / np.linalg.norm(qw[0:4] + q4)
    return tmp[4 - a:4], qw[a:a + 3] + w


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
        # self.w = self.U @ (self.Om - self.w_hkw_vec)
    if len(a_np.shape) == 2:
        return U @ a_np @ U.T
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v]")

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
    raise ValueError(f"Необходимо подать вектор или матрицу! Тип вектора {vec_type} должен быть из [r, v]")


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

        self.to_delete = 0.  # Что ты за попугай?
        self.show_rate = 1  # line_difference в my_plot, line 63

        # Инициализация фильтра
        self.k = KalmanFilter(f=f, c=c, p=self)

        # Инициализация траектории kiam-astro
        if 'kiamastro' in self.v.SOLVER:
            self.jd0 = kiam.juliandate(2024, 1, 1, 0, 0, 0)  # (год, месяц, день, чч, мм, сс)
            self.tr = [[Trajectory(initial_state=np.zeros(6), initial_time=0, initial_jd=self.jd0, variables='rv',
                                   system='gcrs', units_name='earth') for _ in range(obj.n)]
                       for obj in [self.c, self.f, self.a]]
            for j in range(3):
                obj = [self.c, self.f, self.a][j]
                for i in range(obj.n):
                    s0 = np.append(obj.r_irf[i] / (kiam.units('earth')['DistUnit'] * 1e3),
                                   obj.v_irf[i] / (kiam.units('earth')['VelUnit'] * 1e3))
                    self.tr[j][i] = Trajectory(initial_state=s0, initial_time=0, initial_jd=self.jd0, variables='rv',
                                               system='gcrs', units_name='earth')

        # Запись параметров
        self.record = DataFrame()
        self.do_report()  # Продублировать в конце time_step()
        self.record = self.record.astype({'i': 'int32', 'FemtoSat n': 'int32', 'CubeSat n': 'int32'})

    def do_report(self):
        i = self.iter
        d = self.record
        d.loc[i, f'i'] = self.iter
        d.loc[i, f't'] = self.t
        for obj in self.spacecrafts_cd:
            d.loc[i, f'{obj.name} n'] = obj.n
            for i in range(obj.n):
                for v in ['r', 'q', 'v', 'w']:
                    tmp = {'r': [obj.r_irf[i], obj.r_orf[i]], 'v': [obj.v_irf[i], obj.v_orf[i]],
                           'q': [obj.q[i][1:4]], 'w': [obj.w_irf[i], obj.w_orf[i]]}[v]
                    for i_1, frame in enumerate(['orf', 'irf'] if v != 'q' else ['irf']):  # Кватернионы только в ИСК
                        for i_2, c in enumerate('xyz'):
                            d.loc[i, f'{obj.name} {v} {c} {frame} {i}'] = tmp[i_1][i_2]

    # Шаг по времени
    def time_step(self):
        self.iter += 1
        self.t = self.iter * self.v.dT

        # Движение системы
        if 'rk4' in self.v.SOLVER:
            for obj in self.spacecrafts_all:
                for i in range(obj.n):
                    # Вращательное движение
                    obj.q[i], obj.w_irf[i] = rk4_attitude(v_=self.v, obj=obj, i=i)

                    # Поступательное движение
                    if self.v.DYNAMIC_MODEL['aero drag']:
                        obj.r_orf[i], obj.v_orf[i] = rk4_translate(v_=self.v, obj=obj, i=i)
                    else:
                        obj.r_orf[i] = r_hkw(obj.c_hkw[i], self.v.W_ORB, self.t)
                        obj.v_orf[i] = v_hkw(obj.c_hkw[i], self.v.W_ORB, self.t)

                    U, _, _, _ = get_matrices(v=self.v, t=self.t, obj=obj, n=i)
                    obj.r_irf[i] = o_i(v=self.v, a=obj.r_orf[i], U=U, vec_type='r')
                    obj.v_irf[i] = o_i(v=self.v, a=obj.v_orf[i], U=U, vec_type='v')
                    obj.w_orf[i] = i_o(v=self.v, a=obj.w_irf[i], U=U, vec_type='w')
        elif 'kiamastro' in self.v.SOLVER:
            # Расчёт
            if self.iter == 1:
                for i, obj in enumerate(self.spacecrafts_all):
                    for j in range(obj.n):
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

            # Запись положений КА и расчёт вращательного движения КА
            for i, obj in enumerate(self.spacecrafts_all):
                for j in range(obj.n):
                    # Вращательное движение
                    obj.q[j], obj.w_orf[j] = rk4_attitude(v_=self.v, obj=obj, i=j)

                    # Поступательное движение
                    obj.r_irf[j] = np.array([self.tr[i][j].states[ii][self.iter - 1]
                                             for ii in range(3)]) * kiam.units('earth')['DistUnit'] * 1e3
                    obj.v_irf[j] = np.array([self.tr[i][j].states[ii + 3][self.iter - 1]
                                             for ii in range(3)]) * kiam.units('earth')['VelUnit'] * 1e3
                    tr_time = self.tr[i][j].times[self.iter-1] * self.v.SEC_IN_RAD
                    U, _, _, _ = get_matrices(v=self.v, t=tr_time, obj=obj, n=j)
                    obj.r_orf[j] = i_o(v=self.v, a=obj.r_irf[j], U=U, vec_type='r')
                    obj.v_orf[j] = i_o(v=self.v, a=obj.v_irf[j], U=U, vec_type='v')
                    obj.w_orf[j] = i_o(v=self.v, a=obj.w_irf[j], U=U, vec_type='w')
                    self.to_delete = tr_time
        else:
            raise ValueError(f"Поменяй солвер! SOLVER={self.v.SOLVER}, должен быть среди {self.v.SOLVERS}")

        # Комплекс первичной информации
        self.v.MEASURES_VECTOR = []
        self.v.MEASURES_VECTOR_NOTES = []
        measure_antennas_power(c=self.c, f=self.f, v=self.v, noise=np.sqrt(self.v.KALMAN_COEF['r']), produce=True)
        measure_magnetic_field(c=self.c, f=self.f, v=self.v, noise=np.sqrt(self.v.KALMAN_COEF['r']))

        # Изменение режимов работы
        guidance(v=self.v, c=self.c, f=self.f, earth_turn=self.t * self.v.W_ORB / 2 / np.pi)

        # Запись параметров
        self.do_report()
        for obj in self.spacecrafts_all:
            for i in range(obj.n):
                obj.line_orf[i] += [obj.r_orf[i][0], obj.r_orf[i][1], obj.r_orf[i][2]]
                obj.line_irf[i] += [obj.r_irf[i][0], obj.r_irf[i][1], obj.r_irf[i][2]]
        if self.iter % self.show_rate == 0:
            for i_f in range(self.f.n):
                if self.f.operating_mode[i_f] != self.v.OPERATING_MODES[-1]:
                    self.f.line_kalman[i_f] += [self.k.r_orf_estimation[i_f][0],
                                                self.k.r_orf_estimation[i_f][1],
                                                self.k.r_orf_estimation[i_f][2]]
                    self.f.line_difference[i_f] += \
                        [np.array(self.k.r_orf_estimation[i_f][0:3] - np.array(self.f.r_orf[i_f]))]
                    if self.v.NAVIGATION_ANGLES:
                        self.f.attitude_difference[i_f] += [self.k.r_orf_estimation[i_f][3:6]
                                                            - np.array(self.f.q[i_f][1:4])]
                        self.f.spin_difference[i_f] += [self.k.r_orf_estimation[i_f][9:12]
                                                        - np.array(self.f.w_orf[i_f])]
                else:
                    self.f.line_kalman[i_f] += [self.v.NO_LINE_FLAG] * 3
                    self.f.line_difference[i_f] += [self.v.NO_LINE_FLAG * np.ones(3)]
                    if self.v.NAVIGATION_ANGLES:
                        self.f.attitude_difference[i_f] += [self.v.NO_LINE_FLAG * np.ones(3)]
                        self.f.spin_difference[i_f] += [self.v.NO_LINE_FLAG * np.ones(3)]
            for i_c in range(self.c.n):
                for i_f in range(self.f.n):
                    if self.f.operating_mode[i_f] != self.v.OPERATING_MODES[-1]:
                        self.c.kalm_dist[i_c][i_f] += [np.linalg.norm(self.f.r_orf[i_f] -
                                                                      self.k.r_orf_estimation[i_f][0:3])]
                    else:
                        self.c.kalm_dist[i_c][i_f] += [self.v.NO_LINE_FLAG]

        # Навигация чипсатов
        if self.v.IF_NAVIGATION:
            navigate(k=self.k)
