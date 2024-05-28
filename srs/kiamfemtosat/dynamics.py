from datetime import datetime
import kiam_astro
from srs.kiamfemtosat.primary_info import *
from srs.kiamfemtosat.gnc_systems import *
from srs.kiamfemtosat.cosmetic import *
from srs.kiamfemtosat.my_plot import *

# >>>>>>>>>>>> Небольшие функции <<<<<<<<<<<<
def get_atm_params(h: float, h_orb: float) -> tuple:
    """https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html"""
    T = -131.21 + 0.00299 * (h + h_orb)
    p = 2.488 * ((T + 273.1) / 216.6) ** -11.388
    rho = p / (0.2869 * (T + 273.1))
    return rho, T, p

def get_orb_acceleration(r: Union[list, np.ndarray], v: Union[list, np.ndarray], w: float) -> np.ndarray:
    return np.array([-2 * w * v[2],
                     -w ** 2 * r[1],
                     2 * w * v[0] + 3 * w ** 2 * r[2]])

def get_full_acceleration(c_resist: float, square: float, r: Union[list, np.ndarray], v: Union[list, np.ndarray],
                          m: float, h_orb: float, v_orb: float, w: float) -> np.ndarray:
    rho = get_atm_params(h=r[2], h_orb=h_orb)[0]
    force = get_orb_acceleration(r=r, v=v, w=w)
    v_real = v + np.array([v_orb, 0, 0])
    force -= v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / m
    return force

def get_integrate(rvs: Union[list, np.ndarray], dt: float, h_orb: float, v_orb: float, w: float) -> object:
    r_ = np.array(rvs[0:3])
    v_ = np.array(rvs[3:6])
    c_ = rvs[6]
    f = get_full_acceleration(c_resist=1.17, square=c_, m=0.03, r=r_, v=v_, h_orb=h_orb, v_orb=v_orb, w=w)
    v_ += f * dt
    r_ += v_ * dt
    return np.append(np.append(r_, v_), c_).tolist()

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

# >>>>>>>>>>>> Класс динамики кубсатов и чипсатов <<<<<<<<<<<<
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

        # Вывод основных параметров
        if self.iter == 1 and IF_ANY_PRINT:
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

        # Комплекс первичной информации
        measure_antennas_power(c=self.c, f=self.f, noise=np.sqrt(self.r_matrix))
        measure_magnetic_field(c=self.c, f=self.f, noise=np.sqrt(self.r_matrix))

        # Изменение режимов работы
        guidance(c=self.c, f=self.f, earth_turn=self.t * self.w_orb / 2 / np.pi)

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
                        if self.k.orientation:
                            self.f.attitude_difference[i_f] += [self.k.r_orf_estimation[i_f][3:7]
                                                                - np.array(self.f.q[i_f])]
                            self.f.spin_difference[i_f] += [self.k.r_orf_estimation[i_f][10:13]
                                                            - np.array(self.f.w_orf[i_f])]
                    else:
                        self.c.kalm_dist[i_c][i_f] += [NO_LINE_FLAG]
                        self.f.line_kalman[i_f] += [NO_LINE_FLAG] * 3
                        self.f.line_difference[i_f] += [NO_LINE_FLAG * np.ones(3)]
                        if self.k.orientation:
                            self.f.attitude_difference[i_f] += [NO_LINE_FLAG * np.ones(4)]
                            self.f.spin_difference[i_f] += [NO_LINE_FLAG * np.ones(3)]

        # Навигация чипсатов
        if self.k.single_femto_filter:
            for i in range(self.f.n):
                if self.f.operating_mode[i] != OPERATING_MODES[-1]:
                    self.k.calc(i)
                else:
                    self.f.z_difference[i] += [NO_LINE_FLAG]
                if self.c.gain_mode != GAIN_MODES[4]:
                    for j_n in range(self.f.n):
                        for j_t in range(9 if self.k.orientation else 3):  # range(self.k.t)
                            tmp = np.append(np.append(self.f.r_orf[j_n], self.f.v_orf[j_n]), list(self.f.q[j_n][1:4]))
                            tmp = np.zeros(9) if self.f.operating_mode[j_n] != OPERATING_MODES[-1] else tmp
                            self.k.sigmas[j_n * self.k.t + j_t] += [np.sqrt(self.k.p_[j_n][j_t][j_t]) * tmp[j_t]]
        else:
            self.k.calc_all()  # Пока что при любом OPERATING_MODES (если весь рой выпал)

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
