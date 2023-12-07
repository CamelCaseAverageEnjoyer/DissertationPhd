from main_objects import *

def do_tolerance_table(t_integrate: float = 1e4, dt: float = 10., aero: bool = False, repeat: int = 3):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import math

    cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
    kalman_coef = {'q': 1e-15, 'p': [1e-8] * 3, 'r': 1e0}
    tolerance_list = [0.8, 0.9, 1.]  # [0.8, 0.85, 0.90, 0.95, 1.]
    n_f_list = [20, 30, 40]
    counter = 0
    time_begin = datetime.now()
    n_total = len(tolerance_list) * len(tolerance_list) * repeat
    res_mean = [[0. for _ in range(len(n_f_list))] for _ in range(len(tolerance_list))]
    res_std = [[0. for _ in range(len(n_f_list))] for _ in range(len(tolerance_list))]

    for i_t in range(len(tolerance_list)):
        for i_n in range(len(n_f_list)):
            for _ in range(repeat):
                # Отображение процесса
                counter += 1
                per = int(math.ceil(30 * counter / n_total))
                real_per = int(math.ceil(100 * counter / n_total))
                print(f"{real_per}% [{'#' * per + ' ' * (30 - per)}], время: {datetime.now() - time_begin}")

                # Сам процесс
                o = Objects(n_c=1, n_f=n_f_list[i_n], kalman_coef=kalman_coef,
                            model_c='1U', dt=dt, start_navigation='near', if_any_print=False,
                            start_navigation_tolerance=tolerance_list[i_t], method_navigation='kalman_filter rvq all')
                o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas'][3]
                o.c.gain_mode = o.f.gain_mode
                o.p.k.gain_mode = o.f.gain_mode
                o.p.is_aero = aero
                o.p.integrate(t=t_integrate)

                tmp_mean = 0.
                tmp_std = 0.
                for i in range(n_f_list[i_n]):
                    tmp_mean += np.array(o.f.line_difference[i]).mean()
                    tmp_std += np.array(o.f.line_difference[i]).std()
                res_mean[i_t][i_n] += tmp_mean / repeat / n_f_list[i_n]
                res_std[i_t][i_n] += tmp_std / repeat / n_f_list[i_n]

    '''tmp_col = {i: f"Точность {int(tolerance_list[i] * 100)}%" for i in range(len(tolerance_list))}
    tmp_row = {i: f"{n_f_list[i]} аппаратов" for i in range(len(n_f_list))}'''
    tmp_col = {i: f"{int(tolerance_list[i] * 100)}%" for i in range(len(tolerance_list))}
    tmp_row = {i: n_f_list[i] for i in range(len(n_f_list))}
    res_mean_df = pd.DataFrame(res_mean).rename(columns=tmp_col, index=tmp_row)
    res_std_df = pd.DataFrame(res_std).rename(columns=tmp_col, index=tmp_row)
    print(f"Средние средние: \n{res_mean_df}\nСредние отклонения: \n{res_std_df}")
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    sns.heatmap(res_mean_df, annot=True, fmt=".1f", ax=axs[0], cmap=cmap)
    sns.heatmap(res_std_df, annot=True, fmt=".2f", ax=axs[1], cmap=cmap)
    axs[0].set_title("Ошибка определения положения")
    axs[1].set_title("Дисперсия ошибки навигации")
    for i in range(2):
        axs[i].set_xlabel("Точность приорной информации")
        axs[i].set_ylabel("Количество чипсатов")
    plt.show()

def params_search(method_navigation: str = 'kalman_filter rv', t_integrate: float = 1e5, dt: float = 10.,
                  aero: bool = False):
    min_discrepancy = 1e10
    count = 0
    best_params = []
    d_list = [1e-5, 1e-4, 3e-4, 1e-3, 3e-3]
    q_list = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    r_list = [1e-4, 1e-2]
    p1_list = [1e-12, 1e-10, 1e-8]
    p2_list = [1e-12, 1e-10, 1e-8]
    # p3_list = [1e-10, 1e-8]
    # [1e-12, 1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    N = len(d_list) * len(q_list) * len(r_list) * len(p1_list) * len(p2_list)  # * len(p3_list)
    time_begin = datetime.now()
    for q in q_list:
        for r in r_list:
            for p1 in p1_list:
                for p2 in p2_list:
                    # for p3 in p3_list:
                    tmp = 0
                    count += 1
                    if count % 10 == 0:
                        print(f"Поиск коэффициентов фильтра Калмана {count}/{N}={int(100 * count / N)}%, " +
                              real_workload_time(n=count, n_total=N,
                                                 time_begin=time_begin, time_now=datetime.now()))
                    for _ in range(1):
                        o = Objects(n_c=1, n_f=1,
                                    kalman_coef={'q': q, 'p': [p1, p2, p2], 'r': r},
                                    model_c='1U', dt=dt, start_navigation='near',  # random perfect
                                    start_navigation_tolerance=0.98,
                                    method_navigation=method_navigation, if_any_print=False)
                        o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas'][3]
                        o.f.ellipsoidal_signal = 0.99
                        o.p.is_aero = aero
                        o.p.integrate(t=t_integrate)

                        tmp += np.sum(np.abs(np.array(o.c.real_dist[0][0]) - np.array(o.c.calc_dist[0][0])))
                    if tmp < min_discrepancy:
                        min_discrepancy = tmp
                        best_params = [q, p1, p2, p2, r]
    print(f"Подходящие параметры: 'q': {best_params[0]}, 'p': {best_params[1:4]}, 'r': {best_params[4]}")
    print(f"best_params={best_params}")
    return best_params


# params_search()
do_tolerance_table()
