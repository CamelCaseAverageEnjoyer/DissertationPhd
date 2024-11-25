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

def sympy_append(a, b):
    from sympy import Matrix, BlockMatrix
    return Matrix(BlockMatrix([a.T, b.T]).T)

def sympy_mean(a):
    return sum(list(a)) / len(a)

def numerical_and_symbolic_polymorph(trigger_var, trigger_type, trigger_out, not_trigger_out=None):
    def actual_decorator(func):
        def wrapper(*args, **kwargs):
            # print(f"args: {args} | kwargs: {kwargs}")
            trigger = kwargs[trigger_var[1]] if trigger_var[1] in kwargs.keys() else args[trigger_var[0]]
            if isinstance(trigger, trigger_type):
                from numpy import sin, cos, sqrt, tan, arctan as atan, append, pi, mean
                from numpy.linalg import norm, inv
                out_type = trigger_out
                vec_type = np.array
            else:
                from sympy import sin, cos, sqrt, Matrix, atan, tan, pi
                norm = sympy_norm
                append = sympy_append
                mean = sympy_mean
                out_type = Matrix if not_trigger_out is None else not_trigger_out
                vec_type = Matrix
                inv = lambda x: x.inv()
            kwargs['pi'] = pi
            kwargs['sin'] = sin
            kwargs['cos'] = cos
            kwargs['tan'] = tan
            kwargs['sqrt'] = sqrt
            kwargs['norm'] = norm
            kwargs['inv'] = inv
            kwargs['atan'] = atan
            kwargs['append'] = append
            kwargs['mean'] = mean
            kwargs['out_type'] = out_type
            kwargs['vec_type'] = vec_type

            value = func(*args, **kwargs)

            if isinstance(value, tuple):
                return tuple(out_type(i) for i in value)
            return out_type(value)
        return wrapper
    return actual_decorator

'''def sympy_numpy_polymorphism(func):
    """
    Используйте в kwargs аргумент и тип: trigger_var = a, trigger_type = list,
    """
    def local_func(*args, **kwargs):
        if isinstance(kwargs['trigger_var'], kwargs['trigger_type']):
            from numpy import sin, cos, sqrt, tan, arctan as atan
            from numpy.linalg import norm, inv
            kwargs['out'] = kwargs['trigger_out']
        else:
            from sympy import sin, cos, sqrt, Matrix, atan, tan
            norm = sympy_norm
            kwargs['out'] = Matrix
            inv = lambda x: x.inv()
        kwargs['sin'] = sin
        kwargs['cos'] = cos
        kwargs['tan'] = tan
        kwargs['sqrt'] = sqrt
        kwargs['norm'] = norm
        kwargs['inv'] = inv
        kwargs['atan'] = atan

        value = func(*args, **kwargs)

        if isinstance(value, tuple):
            return tuple(kwargs['out'](i) for i in value)
        return kwargs['out'](value)

    return local_func'''

def get_same_type_conversion(a):
    if isinstance(a, list):
        return list
    if isinstance(a, np.ndarray):
        return np.array
    from sympy import Matrix
    if isinstance(a, Matrix):
        return Matrix
