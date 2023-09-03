import numpy as np

def real_workload_time(n: int, n_total: int, time_begin, time_now):
    n_remain = n_total - n
    return f"время: {time_now - time_begin}, оставшееся время: {(time_now - time_begin) * n_remain / n}"

def quart2dcm(L):
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

def clip(a, bot, top):
    if a < bot:
        return bot
    if a > top:
        return top
    return a

def flatten(lst):
    """Функция берёт 2D массив, делает 1D"""
    return [item for sublist in lst for item in sublist]

# РУДИМЕНТЫ
def get_c_hkw(r, v, w):  # Рудимент
    """Возвращает константы C[0]..C[5] движения Хилла-Клохесси-Уилтштира"""
    return [2*r[2] + v[0]/w, v[2]/w, -3*r[2] - 2*v[0]/w, r[0] - 2*v[2]/w, v[1]/w, r[1]]

def r_hkw(C, w, t):  # Рудимент
    """Возвращает вектор координат в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    return np.array([-3*C[0]*w*t + 2*C[1]*np.cos(w*t) - 2*C[2]*np.sin(w*t) + C[3],
                     C[5]*np.cos(w*t) + C[4]*np.sin(w*t),
                     2*C[0] + C[2]*np.cos(w*t) + C[1]*np.sin(w*t)])

def v_hkw(C, w, t):  # Рудимент
    """Возвращает вектор скоростей в момент времени t; \n
    Уравнения движения Хилла-Клохесси-Уилтштира; \n
    Константы C передаются массивом C[0]..C[5]; \n
    Частота w, время t должны быть скалярными величинами."""
    return np.array([-3 * C[0] * w - 2 * w * C[1] * np.sin(w * t) - 2 * w * C[2] * np.cos(w * t),
                     w * C[4] * np.cos(w * t) - w * C[5] * np.sin(w * t),
                     -w * C[2] * np.sin(w * t) + w * C[1] * np.cos(w * t)])
