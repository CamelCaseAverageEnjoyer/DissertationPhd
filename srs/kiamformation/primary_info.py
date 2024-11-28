"""Комплекс первичной информации"""
from spacecrafts import *
from symbolic import numerical_and_symbolic_polymorph


@numerical_and_symbolic_polymorph(trigger_var=(6, 'estimated_params'), trigger_type=(np.ndarray, list),
                                  trigger_out=lambda x: x, not_trigger_out=lambda x: x)
def measure_antennas_power(c: CubeSat, f: FemtoSat, v: Variables, noise: float = None, produce: bool = False,
                           j: int = None, estimated_params=np.array([]), p: any = None, t=None,
                           **kwargs) -> Union[None, tuple]:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры calc_dist при produce==True. Иначе:
    :param c: Класс кубсатов
    :param f: Класс чипсатов
    :param v: Класс гиперпараметров моделирования
    :param noise:
    :param produce: Флаг, указывающий, надо ли записывать полученные величины в PhysicModel.record
    :param j: Количество параметров на 1 дочерний КА
    :param estimated_params: в одну строку
    :param p: Класс PhysicModel (для флага produce)
    :param t: Время (для символьного вычисления)
    :return: None если produce==True (проведение численного моделирования), иначе список измерений + пометки"""
    norm, sqrt, mean, vec_type = kwargs['norm'], kwargs['sqrt'], kwargs['mean'], kwargs['vec_type']
    randy = np.random.uniform(-1, 1, 3)
    anw, notes = [], []
    S_1, S_2, dr, distance = None, None, None, None
    t = p.t if t is None else t

    def get_U(obj, i, t):
        from dynamics import get_matrices
        U, S, A, R_orb = get_matrices(v=v, t=t, obj=obj, n=i)
        return U

    for obj1 in [c, f]:
        for obj2 in [f]:
            for i_1 in range(obj1.n):
                for i_2 in range(obj2.n) if obj1 == c else range(i_1):
                    # >>>>>>>>>>>> Расчёт положений и ориентаций <<<<<<<<<<<<
                    if produce:
                        dr = obj2.r_orf[i_2] - obj1.r_orf[i_1]
                        p.record.loc[p.iter, f'{obj1.name}-{obj2.name} RealDistance {i_1} {i_2}'] = np.linalg.norm(dr)
                        S_1 = quart2dcm(obj1.q[i_1]) @ get_U(obj1, i_1, t).T
                        S_2 = quart2dcm(obj2.q[i_2]) @ get_U(obj2, i_2, t).T
                    else:
                        r1 = vec_type(estimated_params[i_1 * j + 0: i_1 * j + 3]) if obj1 == f else obj1.r_orf[i_1]
                        r2 = vec_type(estimated_params[i_2 * j + 0: i_2 * j + 3])
                        dr = r2 - r1
                        if v.NAVIGATION_ANGLES:
                            q1 = vec2quat(vec_type(estimated_params[i_1 * j + 3: i_1 * j + 6])) \
                                if obj1 == f else obj1.q[i_1]
                            q2 = vec2quat(vec_type(estimated_params[i_2 * j + 3: i_2 * j + 6]))
                            S_1 = quart2dcm(q1) @ get_U(obj1, i_1, t).T
                            S_2 = quart2dcm(q2) @ get_U(obj2, i_2, t).T
                    distance_measured = norm(dr)

                    # >>>>>>>>>>>> Расчёт G и сигнала <<<<<<<<<<<<
                    for direction in ["1->2"]:  # , "2->1"]:
                        take_len = len(get_gain(v=v, obj=obj2 if direction == "1->2" else obj1, r=randy, if_take=True))
                        send_len = len(get_gain(v=v, obj=obj2 if direction == "2->1" else obj1, r=randy, if_send=True))
                        if produce or v.NAVIGATION_ANGLES:
                            G1 = get_gain(v=v, obj=obj1, r=S_1 @ dr,
                                          if_take=direction == "2->1", if_send=direction == "1->2")
                            G2 = get_gain(v=v, obj=obj2, r=S_2 @ dr,
                                          if_take=direction == "1->2", if_send=direction == "2->1")
                            # g_vec = [G1[i] * G2[j] for i in range(take_len if direction == "2->1" else send_len)
                            #                        for j in range(take_len if direction == "1->2" else send_len)]
                            g_vec = [g1 * g2 for g1 in G1 for g2 in G2]
                        else:
                            g_vec = [1] * (take_len if direction == "2->1" else send_len) * \
                                          (take_len if direction == "1->2" else send_len)

                        local_noise = lambda: np.random.normal(0, noise) if produce else 0
                        estimates = [(distance_measured + local_noise()) / sqrt(gg) for gg in g_vec]
                        est_dr = mean(estimates)
                        anw.extend(estimates)

                        # >>>>>>>>>>>> Запись <<<<<<<<<<<<
                        o_fr, i_fr = (obj1, i_1) if direction == "1->2" else (obj2, i_2)
                        o_to, i_to = (obj1, i_1) if direction == "2->1" else (obj2, i_2)
                        if produce:
                            p.record.loc[p.iter, f'{o_fr.name}-{o_to.name} EstimateDistance {i_fr} {i_to}'] = est_dr
                            p.record.loc[p.iter, f'{o_fr.name}-{o_to.name} ErrorEstimateDistance {i_fr} {i_to}'] = \
                                abs(est_dr - norm(dr))
                            p.record.loc[p.iter, f'{o_fr.name}-{o_to.name} ErrorEstimateDistance 1 {i_fr} {i_to}'] = \
                                abs(np.min(estimates) - norm(dr))
                            p.record.loc[p.iter, f'{o_fr.name}-{o_to.name} ErrorEstimateDistance 2 {i_fr} {i_to}'] = \
                                abs(np.max(estimates) - norm(dr))
                        notes.extend([f"{'c' if o_fr == c else 'f'}{'c' if o_to == c else 'f'} {i_fr} {i_to} {j} {i}"
                                      f" {send_len} {take_len}" for i in range(take_len) for j in range(send_len)])

    if produce:
        v.MEASURES_VECTOR = kwargs['vec_type'](anw)
        v.MEASURES_VECTOR_NOTES = notes
    else:
        return kwargs['vec_type'](anw), notes

def measure_magnetic_field(c: CubeSat, f: FemtoSat, v: Variables, noise: float = 0.) -> None:
    """Функция обновляет для объектов CubeSat и FemtoSat параметры b_env"""
    for obj in [c, f]:
        for i in range(obj.n):
            obj.b_env[i] = np.zeros(3) + np.random.normal(0, noise, 3)
    # v.MEASURES_VECTOR += ....????

def measure_gps(f: FemtoSat, noise: float) -> None:
    """Функция обновляет для объектов FemtoSat параметры _не_введено_"""
    pass
