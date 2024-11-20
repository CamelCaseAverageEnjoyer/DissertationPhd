import numpy as np

def get_vars(name: str, n: int, numb: bool = True):
    """Генерит символьные переменные"""
    from sympy import var, Matrix

    s = ""
    axis = ["x", "y", "z"] if n == 3 else [0, "x", "y", "z"]
    for i in range(n):
        s += f"{name}_{i} " if numb else f"{name}_{axis[i]} "

    return Matrix(var(s, real=True))

def get_func(name: str, n: int, numb: bool = True, t=None):
    """Генерит символьные функции"""
    from sympy import Function

    axis = ["x", "y", "z"] if n == 3 else [0, "x", "y", "z"]

    return [Function(f"{name}_{i}" if numb else f"{name}_{axis[i]}", real=True)(t) for i in range(n)]

def sympy_norm(a):
    from sympy import sqrt
    return sqrt(a.dot(a))

def sympy_numpy_polymorphism(func):
    """
    Используйте в kwargs аргумент и тип: trigger_var = a, trigger_type = list,
    """
    def local_func(*args, **kwargs):
        if isinstance(kwargs['trigger_var'], kwargs['trigger_type']):
            from numpy import sin, cos, sqrt
            from numpy.linalg import norm as norm
            kwargs['out'] = kwargs['trigger_out']
        else:
            from sympy import sin, cos, sqrt, Matrix
            norm = sympy_norm
            kwargs['out'] = Matrix
        kwargs['sin'] = sin
        kwargs['cos'] = cos
        kwargs['sqrt'] = sqrt
        kwargs['norm'] = norm

        value = func(*args, **kwargs)

        return kwargs['out'](value)

    return local_func
