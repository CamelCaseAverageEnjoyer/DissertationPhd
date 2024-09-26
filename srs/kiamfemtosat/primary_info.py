"""Комплекс первичной информации. Функции ничего не возвращают, так как запись измерений проводится в элементах
классов КА."""
from srs.kiamfemtosat.spacecrafts import *


def measure_antennas_power(c: CubeSat, f: FemtoSat, v: Variables, noise: float = None, produce: bool = False,
                           get_ready_measurements: bool = False, get_signal_rates: bool = False, 
                           get_model_state: bool = False, j: int = None, x_m: np.ndarray = None) -> Union[None, tuple]:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры calc_dist при produce==True. Иначе:
    1. При get_ready_measurements достаёт значения из v.MEASURES_VECTOR
    2. При get_signal_rates вычисляет G₁G₂
    :param c:
    :param f:
    :param v:
    :param noise:
    :param produce:
    :param get_ready_measurements:
    :param get_signal_rates:
    :param get_model_state: Измерения согласно модели
    :param j:
    :param x_m:
    :return: None если produce==True (проведение численного моделирования), иначе список измерений + пометки"""
    randy = np.random.uniform(-1, 1, 3)
    anw_cf = []
    note_cf = []
    count = 0  # Для get_ready_measurements
    A_1, A_2, A_f, A_c, dr, dist_estimate = None, None, None, None, None, None
    for i_c in range(c.n):
        for i_f in range(f.n):
            if produce:
                dr = f.r_orf[i_f] - c.r_orf[i_c]
                c.real_dist[i_c][i_f] += [np.linalg.norm(dr)]
                A_f = quart2dcm(np.array(f.q[i_f]))
                A_c = quart2dcm(np.array(c.q[i_c]))
                dist_estimate = np.random.normal(0, noise) + np.linalg.norm(dr)
            elif get_signal_rates or get_model_state:
                dr = x_m[i_f*j + 0: i_f*j + 3] - c.r_orf[i_c]
                if get_signal_rates:
                    A_f = quart2dcm(x_m[i_f*j + 3: i_f*j + 7])
                    A_c = quart2dcm(c.q[i_c])

            for direction in ["f->c", "c->f"]:
                take_len = len(get_gain(v=v, obj=f if direction == "c->f" else c, r=randy, if_take=True))
                send_len = len(get_gain(v=v, obj=f if direction == "f->c" else c, r=randy, if_send=True))
                if produce or get_signal_rates:
                    anw_tmp = flatten([[
                        get_gain(v=v, obj=f, r=A_f @ dr, if_take=direction == "c->f", if_send=direction == "f->c")[ii] *
                        get_gain(v=v, obj=c, r=A_c @ dr, if_take=direction == "f->c", if_send=direction == "c->f")[jj]
                        for ii in range(take_len if direction == "c->f" else send_len)]
                        for jj in range(take_len if direction == "f->c" else send_len)])
                    if produce:
                        anw_tmp = [dist_estimate/np.sqrt(i_a) for i_a in anw_tmp]
                elif get_ready_measurements:
                    anw_tmp = v.MEASURES_VECTOR[count: count + take_len*send_len]
                    count += take_len*send_len
                elif get_model_state:
                    anw_tmp = [np.linalg.norm(dr) for _ in range(take_len * send_len)]
                else:
                    anw_tmp = None

                if direction == "f->c":
                    note_cf += flatten([[f"fc {i_f} {i_c} {ii} {jj} {send_len} {take_len}"
                                         for ii in range(send_len)] for jj in range(take_len)])
                else:
                    c.dist_estimate[i_c][i_f] += [np.mean(anw_tmp)] if produce else []
                    note_cf += flatten([[f"cf {i_c} {i_f} {jj} {ii} {send_len} {take_len}"
                                         for ii in range(take_len)] for jj in range(send_len)])
                anw_cf += anw_tmp

    anw_ff = []
    note_ff = []
    for i_f1 in range(f.n):
        for i_f2 in range(i_f1):
            if produce:
                dr = f.r_orf[i_f1] - f.r_orf[i_f2]
                f.real_dist[i_f1][i_f2] += [np.linalg.norm(dr)]
                A_1 = quart2dcm(np.array(f.q[i_f1]))
                A_2 = quart2dcm(np.array(f.q[i_f2]))
                dist_estimate = np.random.normal(0, noise) + np.linalg.norm(dr)
            elif get_signal_rates or get_model_state:
                dr = x_m[i_f1*j + 0: i_f1*j + 3] - x_m[i_f2*j + 0: i_f2*j + 3]
                if get_signal_rates:
                    A_1 = quart2dcm(x_m[i_f1*j + 3: i_f1*j + 6])
                    A_2 = quart2dcm(x_m[i_f2*j + 3: i_f2*j + 6])

            for direction in ["2->1", "1->2"]:
                take_len = len(get_gain(v=v, obj=f, r=randy, if_take=True))
                send_len = len(get_gain(v=v, obj=f, r=randy, if_send=True))
                if produce or get_signal_rates:
                    anw_tmp = flatten([[
                        get_gain(v=v, obj=f, r=A_1 @ dr, if_take=direction == "2->1", if_send=direction == "1->2")[ii] *
                        get_gain(v=v, obj=f, r=A_2 @ dr, if_take=direction == "1->2", if_send=direction == "2->1")[jj]
                        for ii in range(take_len if direction == "2->1" else send_len)]
                        for jj in range(take_len if direction == "1->2" else send_len)])
                    if produce:
                        anw_tmp = [dist_estimate/np.sqrt(i_a) for i_a in anw_tmp]
                elif get_ready_measurements:
                    anw_tmp = v.MEASURES_VECTOR[count: count + take_len*send_len]
                    count += take_len*send_len
                elif get_model_state:
                    anw_tmp = [np.linalg.norm(dr) for _ in range(take_len * send_len)]
                else:
                    anw_tmp = None

                anw_ff += anw_tmp
                if direction == "1->2":
                    f.dist_estimate[i_f1][i_f2] += [np.mean(anw_tmp)] if produce else []
                    note_ff += flatten([[f"ff {i_f1} {i_f2} {jj} {ii} {send_len} {take_len}"
                                         for ii in range(take_len)] for jj in range(send_len)])
                else:
                    f.dist_estimate[i_f2][i_f1] += [np.mean(anw_tmp)] if produce else []
                    note_ff += flatten([[f"ff {i_f2} {i_f1} {jj} {ii} {send_len} {take_len}"
                                         for ii in range(take_len)] for jj in range(send_len)])

    if produce:
        v.MEASURES_VECTOR += anw_cf
        # my_print(f"Длина измерений 1: {len(v.MEASURES_VECTOR)}", color='r', if_print=v.IF_TEST_PRINT)
        v.MEASURES_VECTOR += anw_ff
        # my_print(f"Длина измерений 2: {len(v.MEASURES_VECTOR)}", color='r', if_print=v.IF_TEST_PRINT)
        v.MEASURES_VECTOR_NOTES += note_cf
        v.MEASURES_VECTOR_NOTES += note_ff
    else:
        return anw_cf, anw_ff, note_cf + note_ff

def measure_magnetic_field(c: CubeSat, f: FemtoSat, v: Variables, noise: float = 0.) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            obj.b_env[i] = np.zeros(3) + np.random.normal(0, noise, 3)
    # v.MEASURES_VECTOR += ....????

def measure_gps(f: FemtoSat, noise: float) -> None:
    """Функция обновляет для объектов FemtoSat параметры _не_введено_"""
    pass
