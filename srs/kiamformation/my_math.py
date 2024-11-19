"""Небольшие функции"""
from typing import Union
import numpy as np
import quaternion


def get_antisymmetric_matrix(a):
    return np.array([[0, -a[2], a[1]],
                     [a[2], 0, -a[0]],
                     [-a[1], a[0], 0]])

def deg2rad(a: float) -> float:
    return a / 180 * np.pi

def rad2deg(a: float) -> float:
    return a * 180 / np.pi

def vec2unit(a: np.ndarray) -> np.ndarray:
    return a / np.linalg.norm(a)

def my_cross(a, b):
    """Функция векторного произведения"""
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])

def vec2quat(v):
    """Перевод вектора кватерниона в кватернион"""
    if isinstance(v, np.ndarray):
        return np.quaternion(np.sqrt(1 - np.linalg.norm(v)**2), *v)
    else:
        from sympy import Matrix, sqrt
        return Matrix([sqrt(1 - v.dot(v)), v[0], v[1], v[2]])

def q_dot(q1, q2):
    """Умножение кватернионов длины 4"""
    if isinstance(q1, quaternion.quaternion):
        return q1 * q2
    else:
        from sympy import Matrix
        q1v = Matrix(q1[1:4])
        q2v = Matrix(q2[1:4])
        scalar = q1[0]*q2[0] - q1v.dot(q2v)
        vector = q1v*q2[0] + q2v*q1[0] + q1v.cross(q2v)
        return Matrix([scalar, vector[0], vector[1], vector[2]])

def euler2rot_matrix(a: float, b: float, g: float) -> np.ndarray:
    return quaternion.as_rotation_matrix(quaternion.from_euler_angles(a, b, g))  # А где ты вообще нужен?

def pol2dec(r, u, v):
    return np.array([r * np.cos(u) * np.cos(v), r * np.sin(u) * np.cos(v), r * np.sin(v)])

def quart2dcm(q):
    """Матрицу поворота из кватерниона поворота"""
    if isinstance(q, quaternion.quaternion):
        return quaternion.as_rotation_matrix(q)
    else:
        from sympy import Matrix
        w, x, y, z = q
        return Matrix([[1 - 2*y**2 - 2*z**2, 2*x*y + 2*z*w, 2*x*z - 2*y*w],
                       [2*x*y - 2*z*w, 1 - 2*x**2 - 2*z**2, 2*y*z + 2*x*w],
                       [2*x*z + 2*y*w, 2*y*z - 2*x*w, 1 - 2*x**2 - 2*y**2]])

def clip(a: float, bot: float, top: float) -> float:
    if a < bot:
        return bot
    if a > top:
        return top
    return a

def flatten(lst: Union[list, np.ndarray]) -> list:
    """Функция берёт 2D массив, делает 1D"""
    return [item for sublist in lst for item in sublist]
