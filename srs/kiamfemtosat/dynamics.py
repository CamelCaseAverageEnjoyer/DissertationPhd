"""Функции для моделирования динамики КА"""
from datetime import datetime
from kiam_astro import kiam
from kiam_astro.trajectory import Trajectory

from srs.kiamfemtosat.cosmetic import *
from srs.kiamfemtosat.primary_info import *
from srs.kiamfemtosat.gnc_systems import *
from srs.kiamfemtosat.my_plot import *

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
def get_atm_params(h: float, atm_model: str = None) -> tuple:
    """NASA модель: https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html (якорная точка: 25 км)
    ПНБО модель: https://www.energia.ru/ktt/archive/2022/04-2022/101-111.pdf (120-600 км)
    COESA62, COESA76 модели: библиотека poliastro
    :param h: Высота
    :param atm_model: Модель атмосферы (необязательный параметр)
    :return: (ρ, T, P): плотность, температура, давление (Внимание! Для ПНБО (ρ, None, None)"""
    atm_model = ATMOSPHERE_MODEL if atm_model is None else atm_model
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

def get_geopotential_acceleration(r: Union[list, np.ndarray], v: Union[list, np.ndarray]) -> np.ndarray:
    """Возвращает ускорение КА от притяжения Земли.
    Внимание! При глобальном параметре DYNAMIC_MODEL='Clohessy-Wiltshire' возвращает ускорение в ОСК.
    Иначе возвращает ускорение в ИСК."""
    if DYNAMIC_MODEL == 'Clohessy-Wiltshire':
        return np.array([-2 * W_ORB * v[2],
                         -W_ORB ** 2 * r[1],
                         2 * W_ORB * v[0] + 3 * W_ORB ** 2 * r[2]])
    return MU * np.array(r) / np.linalg.norm(r)**3

def get_aero_drag_acceleration(obj: Union[CubeSat, FemtoSat], i: int, r: Union[list, np.ndarray],
                               v: Union[list, np.ndarray]):
    """Возвращает ускорение КА от сопротивления атмосферы.
    Внимание! При глобальном параметре DYNAMIC_MODEL='Clohessy-Wiltshire' возвращает ускорение в ОСК.
    Иначе возвращает ускорение в ИСК."""
    S = quart2dcm(obj.q[i])  # Для новых кватернионов - неверно!
    cos_alpha = clip((np.trace(S) - 1) / 2, -1, 1)
    # alpha = 180 / np.pi * np.arccos(cos_alpha)
    rho = get_atm_params(h=obj.r_orf[i][2] + ORBIT_RADIUS - EARTH_RADIUS)[0]
    if obj.name == "CubeSat":
        c_resist = 1.05
        square = obj.size[0] * obj.size[1]
    else:
        c_resist = 1.17
        square = obj.size[0] * obj.size[1] * abs(cos_alpha)

    if DYNAMIC_MODEL == 'Clohessy-Wiltshire':
        v_real = v + np.array([V_ORB, 0, 0])
        rho = get_atm_params(h=r[2] + ORBIT_RADIUS - EARTH_RADIUS)[0]
        return - v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / obj.mass

def get_full_acceleration(obj: Union[CubeSat, FemtoSat], i: int, r: Union[float, np.ndarray],
                          v: Union[float, np.ndarray]) -> np.ndarray:
    """Возвращает вектор силы в ОСК, принимает параметры в ОСК\n
    square - площадь S для аэродинамики
    c_resist - Cf для аэродинаммики"""
    if DYNAMIC_MODEL == 'Clohessy-Wiltshire':
        force = get_geopotential_acceleration(r, v)
        if AERO_DRAG:
            force += get_aero_drag_acceleration(r=r, v=v, obj=obj, i=i)
        return force

def rk4_translate(obj: Union[CubeSat, FemtoSat], i: int, dt: float = dT, r=None, v=None) -> tuple:
    def rv_right_part(rv1, a1):
        return np.array([rv1[3], rv1[4], rv1[5], a1[0], a1[1], a1[2]])
    r = obj.r_orf[i] if r is None else r
    v = obj.v_orf[i] if v is None else v
    a = get_full_acceleration(obj=obj, i=i, r=r, v=v)

    rv = np.append(r, v)
    k1 = rv_right_part(rv, a)
    k2 = rv_right_part(rv + k1 * dt / 2, a)
    k3 = rv_right_part(rv + k2 * dt / 2, a)
    k4 = rv_right_part(rv + k3 * dt, a)
    rv = dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return rv[0:3] + r, rv[3:6] + v


# >>>>>>>>>>>> Вращательное движение, интегрирование <<<<<<<<<<<<
def get_torque(obj: Union[CubeSat, FemtoSat], q: Union[list, np.ndarray], w: np.ndarray) -> np.ndarray:
    """Возвращает вектор углового УСКОРЕНИЯ. Удачи!"""
    return np.zeros(3)

def rk4_attitude(obj: Union[CubeSat, FemtoSat], i: int, dt: float = dT, q=None, w=None):
    """Господи, где здесь производная, где дифференциал? Пожалуйста, дрогой я, дай мне знак!"""
    def lw_right_part(qw_, e_):
        q_, w_ = qw_[0:4], qw_[4:7]
        dq = 1 / 2 * q_dot([0, w_[0], w_[1], w_[2]], q_)
        return np.append(dq, e_)
    q = obj.q[i] if q is None else q
    w = obj.w_orf[i] if w is None else w
    e = get_torque(obj=obj, q=q, w=w)

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

# >>>>>>>>>>>> Класс динамики кубсатов и чипсатов <<<<<<<<<<<<
class PhysicModel:
    def __init__(self, f: FemtoSat, c: CubeSat):
        self.show_rate = 1

        # Неизменные параметры
        self.t = 0.
        self.iter = 0
        self.c = c
        self.f = f
        self.time_begin = datetime.now()

        # Инициализация фильтра
        self.k = KalmanFilter(f=f, c=c, p=self)

        # Расчёт движения Хилла-Клохесси-Уилтшира
        self.c.c_hkw = [get_c_hkw(self.c.r_orf[i], self.c.v_orf[i], W_ORB) for i in range(self.c.n)]
        self.f.c_hkw = [get_c_hkw(self.f.r_orf[i], self.f.v_orf[i], W_ORB) for i in range(self.f.n)]

    # Шаг по времени
    def time_step(self):
        self.iter += 1
        self.t = self.iter * dT

        if self.iter == 1 and IF_ANY_PRINT:
            # Отображение траектории полёта
            t0 = 0.0
            s0 = [ORBIT_RADIUS / EARTH_RADIUS, ECCENTRICITY, kiam.deg2rad(INCLINATION),
                  0.0, 0.0, 0.0]  # a, e, i, Ω, ω, M₀
            # s0.extend(list(kiam.eye2vec(6)))  # Orbital_elements + State-transition matrix
            s0 = np.array(s0)
            jd0 = kiam.juliandate(2024, 1, 1, 0, 0, 0)  # (год, месяц, день, чч, мм, сс)
            tr = Trajectory(initial_state=s0, initial_time=t0, initial_jd=jd0, variables='oe', system='gcrs',
                            units_name='earth')
            tr.set_model(variables='rv', model_type='nbp', primary='earth',
                         sources_list=['atm'])  # 'atm', 'j2', 'moon', 'sun'
            tr.model['data']['jd_zero'] = jd0  # julian date corresponding to t = 0
            tr.model['data']['mass'] = self.f.mass  # spacecraft mass, kg
            tr.model['data']['area'] = self.f.size[0] * self.f.size[1]  # spacecraft area, m^2
            tr.model['data']['order'] = 0  # order of the Moon's gravity field
            tr.propagate(tof=2 * np.pi * TIME*W_ORB/2/np.pi, npoints=50000)  # (time of flight, number of points)
            print(f"Time: {tr.times[-1]} ({2 * np.pi})")
            # help(tr.show)
            tmp = tr.show(variables='3d', language='rus')

            # Вывод основных параметров
            tmp = ", ориентации" if NAVIGATION_ANGLES else ""
            my_print(f"Диаграмма антенн кубсата: {self.c.gain_mode}\n"
                     f"Диаграмма антенн фемтосатов: {self.f.gain_mode}\n"
                     f"Учёт аэродинамики: {AERO_DRAG}\n"
                     f"Применяется фильтр Калмана для поправки: положений, скоростей{tmp}\n" 
                     f"Фильтр Калмана основан на: "
                     f"{'всех чипсатах' if NAVIGATION_BY_ALL else 'одном чипсате'}", color='c')
            if not IF_NAVIGATION:
                my_print(f"Внимание: IF_NAVIGATION={IF_NAVIGATION}! ", color='m')

        # Движение системы на dt
        for obj in [self.c, self.f]:
            for i in range(obj.n):
                # Вращательное движение
                obj.q[i], obj.w_orf[i] = rk4_attitude(obj=obj, i=i)

                # Поступательное движение
                if AERO_DRAG:
                    obj.r_orf[i], obj.v_orf[i] = rk4_translate(obj=obj, i=i)
                else:
                    obj.r_orf[i] = r_hkw(obj.c_hkw[i], W_ORB, self.t)
                    obj.v_orf[i] = v_hkw(obj.c_hkw[i], W_ORB, self.t)

                if self.iter % self.show_rate == 0:
                    obj.line[i] += [obj.r_orf[i][0], obj.r_orf[i][1], obj.r_orf[i][2]]

        # Комплекс первичной информации
        measure_antennas_power(c=self.c, f=self.f, noise=np.sqrt(KALMAN_COEF['r']))
        measure_magnetic_field(c=self.c, f=self.f, noise=np.sqrt(KALMAN_COEF['r']))

        # Изменение режимов работы
        guidance(c=self.c, f=self.f, earth_turn=self.t * W_ORB / 2 / np.pi)

        # Запись параметров
        if self.iter % self.show_rate == 0:
            for i_c in range(self.c.n):
                for i_f in range(self.f.n):
                    if self.f.operating_mode[i_f] != OPERATING_MODES[-1]:
                        self.c.kalm_dist[i_c][i_f] += [np.linalg.norm(self.f.r_orf[i_f] -
                                                                      self.k.r_orf_estimation[i_f][0:3])]
                        self.f.line_kalman[i_f] += [self.k.r_orf_estimation[i_f][0],
                                                    self.k.r_orf_estimation[i_f][1],
                                                    self.k.r_orf_estimation[i_f][2]]
                        self.f.line_difference[i_f] += \
                            [np.array(self.k.r_orf_estimation[i_f][0:3] - np.array(self.f.r_orf[i_f]))]
                        if NAVIGATION_ANGLES:
                            self.f.attitude_difference[i_f] += [self.k.r_orf_estimation[i_f][3:7]
                                                                - np.array(self.f.q[i_f])]
                            self.f.spin_difference[i_f] += [self.k.r_orf_estimation[i_f][10:13]
                                                            - np.array(self.f.w_orf[i_f])]
                    else:
                        self.c.kalm_dist[i_c][i_f] += [NO_LINE_FLAG]
                        self.f.line_kalman[i_f] += [NO_LINE_FLAG] * 3
                        self.f.line_difference[i_f] += [NO_LINE_FLAG * np.ones(3)]
                        if NAVIGATION_ANGLES:
                            self.f.attitude_difference[i_f] += [NO_LINE_FLAG * np.ones(4)]
                            self.f.spin_difference[i_f] += [NO_LINE_FLAG * np.ones(3)]

        # Навигация чипсатов
        if IF_NAVIGATION:
            navigate(k=self.k)

    # Перевод между системами координат
    def get_matrices(self, obj: Union[CubeSat, FemtoSat], n: int, t: float = None):
        t = self.t if t is None else t
        A = quart2dcm(obj.q[n])
        U = np.array([[0., 1., 0.],
                      [0., 0., 1.],
                      [1., 0., 0.]]) @ \
            np.array([[np.cos(t * W_ORB), np.sin(t * W_ORB), 0],
                      [-np.sin(t * W_ORB), np.cos(t * W_ORB), 0],
                      [0, 0, 1]])
        S = A @ U.T
        R_orb = U.T @ np.array([0, 0, ORBIT_RADIUS])
        return U, S, A, R_orb
