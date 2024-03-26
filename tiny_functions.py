"""Функции размеров не столь отдалённых"""
import numpy as np
from typing import Union


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

def get_gain(r: Union[float, np.ndarray]) -> float:
    r1 = r / np.linalg.norm(r)
    return 2 \
        - np.linalg.norm(np.dot(r1, np.array([1, 0, 0]))) ** 2 \
        - np.linalg.norm(np.dot(r1, np.array([0, 1, 0]))) ** 2

def get_rand_c(w: float, r_spread: float = 100, v_spread: float = 0.01, if_no: bool = False,
               if_quaternion: bool = False) -> list:
    """(quaternion or quaternion_but_i_dont_give_a_fuck)"""
    x, y, z = np.random.uniform(-r_spread, r_spread, 3)
    vx, vy, vz = np.random.uniform(-v_spread, v_spread, 3)
    a, b, g = np.random.uniform(-100, 100, 3)
    if if_quaternion and not if_no:
        return [2 * z + vx / w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y, a, b, g]
    return [2 * z + vx / w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y]

def get_quaternion_from_vecpart(q: Union[list, np.ndarray]) -> np.ndarray:
    return np.append(1 - np.linalg.norm(q) ** 2, q) if np.linalg.norm(q) <= 1 else \
        np.append(0, q / np.linalg.norm(q))

def real_workload_time(n: int, n_total: int, time_begin, time_now) -> str:
    n_remain = n_total - n
    return f"время: {time_now - time_begin}, оставшееся время: {(time_now - time_begin) * n_remain / n}"

def vec2quat(v: Union[list, np.ndarray]) -> list:
    """Перевод вектор-часть кватерниона в кватернион"""
    if len(v) != 3:
        raise ValueError(f"Подаётся вектор длинны {len(v)}, требуется длинна 3!")
    if np.linalg.norm(v) > 1:
        return [0] + list(np.array(v)/np.linalg.norm(v))
    else:
        return [np.sqrt(1 - np.linalg.norm(v)**2), v[0], v[1], v[2]]

def q_dot(L1, L2):
    """Функция является кватернионным умножением; \n
    Кватернион L1,L2 передаются векторами длины 4; \n
    Возвращает кватернион L[0]..L[3]."""
    return np.array([L1[0] * L2[0] - L1[1] * L2[1] - L1[2] * L2[2] - L1[3] * L2[3],
                     L1[0] * L2[1] + L1[1] * L2[0] + L1[2] * L2[3] - L1[3] * L2[2],
                     L1[0] * L2[2] + L1[2] * L2[0] + L1[3] * L2[1] - L1[1] * L2[3],
                     L1[0] * L2[3] + L1[3] * L2[0] + L1[1] * L2[2] - L1[2] * L2[1]])

def euler2rot_matrix(a: float, b: float, g: float) -> np.ndarray:
    return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]]) @ \
        np.array([[1, 0, 0], [0, np.cos(b), -np.sin(b)], [0, np.sin(b), np.cos(b)]]) @ \
        np.array([[np.cos(g), -np.sin(g), 0], [np.sin(g), np.cos(g), 0], [0, 0, 1]])

def quart2dcm(L) -> np.ndarray:
    """Функция ищет матрицу поворота из кватерниона поворота; \n
    Кватернион L передаётся вектором длины 4; \n
    Возвращает матрицу 3х3."""
    w, x, y, z = L
    A = np.eye(3)
    A[0][0] = 1 - 2 * y ** 2 - 2 * z ** 2
    A[0][1] = 2 * x * y + 2 * z * w
    A[0][2] = 2 * x * z - 2 * y * w
    A[1][0] = 2 * x * y - 2 * z * w
    A[1][1] = 1 - 2 * x ** 2 - 2 * z ** 2
    A[1][2] = 2 * y * z + 2 * x * w
    A[2][0] = 2 * x * z + 2 * y * w
    A[2][1] = 2 * y * z - 2 * x * w
    A[2][2] = 1 - 2 * x ** 2 - 2 * y ** 2
    return A

def clip(a: float, bot: float, top: float) -> float:
    if a < bot:
        return bot
    if a > top:
        return top
    return a

def flatten(lst: Union[list, np.ndarray]) -> list:
    """Функция берёт 2D массив, делает 1D"""
    return [item for sublist in lst for item in sublist]

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
