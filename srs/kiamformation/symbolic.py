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
