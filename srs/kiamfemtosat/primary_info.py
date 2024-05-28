"""Комплекс первичной информации"""
from srs.kiamfemtosat.spacecrafts import *


def measure_antennas_power(c: CubeSat, f: FemtoSat, noise: float) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры calc_dist"""
    for i_c in range(c.n):
        for i_f in range(f.n):
            dr = f.r_orf[i_f] - c.r_orf[i_c]
            c.real_dist[i_c][i_f] += [np.linalg.norm(dr)]
            A_f = quart2dcm(np.array(f.q[i_f]))
            A_c = quart2dcm(np.array(c.q[i_c]))
            calc_dist = np.random.normal(0, noise) + np.linalg.norm(dr) / \
                np.sqrt(get_gain(f, r=A_f @ dr, mode3=True)[0] * get_gain(c, r=A_c @ dr, mode3=True)[0])

            c.calc_dist[i_c][i_f] += [calc_dist]

    for i_f1 in range(f.n):
        for i_f2 in range(f.n):
            if i_f1 != i_f2:
                dr = f.r_orf[i_f1] - f.r_orf[i_f2]
                f.real_dist[i_f1][i_f2] += [np.linalg.norm(dr)]

                q1 = f.q[i_f1]
                q2 = f.q[i_f1]
                calc_dist = np.random.normal(0, noise) + np.linalg.norm(dr) / \
                    np.sqrt(get_gain(f, r=quart2dcm(q1)@dr, mode3=True)[0] *
                            get_gain(f, r=quart2dcm(q2)@dr, mode3=True)[0])
                f.calc_dist[i_f1][i_f2] += [calc_dist]
            else:
                f.calc_dist[i_f1][i_f2] += [0.]

def measure_magnetic_field(c: CubeSat, f: FemtoSat, noise: float = 0.) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            obj.b_env[i] = np.zeros(3) + np.random.normal(0, noise, 3)
