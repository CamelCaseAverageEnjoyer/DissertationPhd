"""Комплекс первичной информации. Функции ничего не возвращают, так как запись измерений проводится в элементах
классов КА."""
from srs.kiamfemtosat.spacecrafts import *


def measure_antennas_power(c: CubeSat, f: FemtoSat, v: Variables, noise: float) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры calc_dist"""
    anw_cf = []
    for i_c in range(c.n):
        for i_f in range(f.n):
            dr = f.r_orf[i_f] - c.r_orf[i_c]
            c.real_dist[i_c][i_f] += [np.linalg.norm(dr)]
            A_f = quart2dcm(np.array(f.q[i_f]))
            A_c = quart2dcm(np.array(c.q[i_c]))
            dist_estimate = np.random.normal(0, noise) + np.linalg.norm(dr)
            # Измерения чипсата сигналов от кубсата
            anw_tmp = []
            take_len = len(get_gain(v=v, obj=f, r=A_f @ dr, if_take=True))
            send_len = len(get_gain(v=v, obj=c, r=A_c @ dr, if_send=True))
            anw_tmp += flatten([[dist_estimate/np.sqrt(get_gain(v=v, obj=f, r=A_f @ dr, if_take=True)[ii] *
                                                       get_gain(v=v, obj=c, r=A_c @ dr, if_send=True)[jj])
                                  for ii in range(take_len)] for jj in range(send_len)])
            v.MEASURES_VECTOR_NOTES += flatten([[f"cf {i_c} {i_f} {jj} {ii} {send_len} {take_len}"
                                                 for ii in range(take_len)] for jj in range(send_len)])
            # Измерения кубсата сигналов от чипсата
            take_len = len(get_gain(v=v, obj=c, r=A_f @ dr, if_take=True))
            send_len = len(get_gain(v=v, obj=f, r=A_c @ dr, if_send=True))
            anw_tmp += flatten([[dist_estimate/np.sqrt(get_gain(v=v, obj=f, r=A_f @ dr, if_send=True)[ii] *
                                                       get_gain(v=v, obj=c, r=A_c @ dr, if_take=True)[jj])
                                   for ii in range(send_len)] for jj in range(take_len)])
            v.MEASURES_VECTOR_NOTES += flatten([[f"fc {i_f} {i_c} {ii} {jj} {send_len} {take_len}"
                                                 for ii in range(send_len)] for jj in range(take_len)])
            anw_cf += anw_tmp
            c.dist_estimate[i_c][i_f] += [np.mean(anw_tmp)]

    anw_ff = []
    for i_f1 in range(f.n):
        for i_f2 in range(f.n):
            if i_f1 != i_f2:
                dr = f.r_orf[i_f1] - f.r_orf[i_f2]
                f.real_dist[i_f1][i_f2] += [np.linalg.norm(dr)]
                A_1 = quart2dcm(np.array(f.q[i_f1]))
                A_2 = quart2dcm(np.array(f.q[i_f2]))
                dist_estimate = np.random.normal(0, noise) + np.linalg.norm(dr)
                # Измерения чипсата 1 сигналов от чипсата 2
                anw_tmp = []
                take_len = len(get_gain(v=v, obj=f, r=A_1 @ dr, if_take=True))
                send_len = len(get_gain(v=v, obj=f, r=A_2 @ dr, if_send=True))
                anw_tmp += flatten([[dist_estimate/np.sqrt(get_gain(v=v, obj=f, r=A_1 @ dr, if_take=True)[ii] *
                                                           get_gain(v=v, obj=f, r=A_2 @ dr, if_send=True)[jj])
                                     for ii in range(take_len)] for jj in range(send_len)])
                v.MEASURES_VECTOR_NOTES += flatten([[f"ff {i_f2} {i_f1} {jj} {ii} {send_len} {take_len}"
                                                     for ii in range(take_len)] for jj in range(send_len)])
                anw_ff += anw_tmp
                f.dist_estimate[i_f1][i_f2] += [np.mean(anw_tmp)]

    v.MEASURES_VECTOR += anw_cf
    # my_print(f"Длина измерений 1: {len(v.MEASURES_VECTOR)}", color='r', if_print=v.IF_TEST_PRINT)
    v.MEASURES_VECTOR += anw_ff
    # my_print(f"Длина измерений 2: {len(v.MEASURES_VECTOR)}", color='r', if_print=v.IF_TEST_PRINT)

def measure_magnetic_field(c: CubeSat, f: FemtoSat, v: Variables, noise: float = 0.) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            obj.b_env[i] = np.zeros(3) + np.random.normal(0, noise, 3)
    # v.MEASURES_VECTOR += ....????

def measure_gps(f: FemtoSat, noise: float) -> None:
    """Функция обновляет для объектов FemtoSat параметры _не_введено_"""
    pass
