import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import rcParams
from srs.kiamfemtosat.my_math import *
from srs.kiamfemtosat.config import *


FEMTO_RATE = 1e2
CUBE_RATE = 1e2
TITLE_SIZE = 15  # 15
CAPTION_SIZE = 13  # 13
rcParams["savefig.directory"] = "/home/kodiak/Desktop"
rcParams["savefig.format"] = "jpg"

def show_chipsat(o, j, clr, opacity):
    global FEMTO_RATE
    x, y, z = ([], [], [0 for _ in range(4)])
    for x_shift in [-o.f.size[0] * FEMTO_RATE, o.f.size[0] * FEMTO_RATE]:
        for y_shift in [-o.f.size[1] * FEMTO_RATE, o.f.size[1] * FEMTO_RATE]:
            x += [x_shift]
            y += [y_shift]
    A = quart2dcm(o.f.q[j])
    for i in range(4):
        r = A.T @ np.array([x[i], y[i], z[i]])
        x[i] = r[0] + o.f.r_orf[j][0]
        y[i] = r[1] + o.f.r_orf[j][1]
        z[i] = r[2] + o.f.r_orf[j][2]
    xl, yl, zl = ([], [], [])
    for i in range(int(len(o.f.line[j])//3)):
        xl += [o.f.line[j][i * 3 + 0]]
        yl += [o.f.line[j][i * 3 + 1]]
        zl += [o.f.line[j][i * 3 + 2]]
    xk, yk, zk = ([], [], [])
    for i in range(int(len(o.f.line_kalman[j])//3)):
        xk += [o.f.line_kalman[j][i * 3 + 0]]
        yk += [o.f.line_kalman[j][i * 3 + 1]]
        zk += [o.f.line_kalman[j][i * 3 + 2]]
    return [go.Mesh3d(x=x, y=y, z=z, color=clr, opacity=opacity), go.Scatter3d(x=xl, y=yl, z=zl, mode='lines'),
            go.Scatter3d(x=xk, y=yk, z=zk, mode='lines')]

def show_cubesat(o, j):
    global CUBE_RATE
    total_cubes = 30
    r = [[[] for _ in range(total_cubes)] for _ in range(3)]
    # r[x/y/z][0..5 - yellow, 6..29 - gray][1..4 - sides of square]
    sequence = [[0, 0], [0, 1], [1, 0], [1, 1]]
    for i in range(3):
        for k in range(2):
            r[i][k + i * 2] = [(-1)**(k+1) * o.c.size[i] * CUBE_RATE for _ in range(4)]

    shift = [[-o.c.size[i] * CUBE_RATE, o.c.size[i] * CUBE_RATE] for i in range(3)]
    legs = [o.c.legs_x, o.c.legs_x, o.c.legs_z]
    n_legs = int((total_cubes - 6) // 6)
    bound_legs = [[[((-1)**sequence[s][0] * o.c.size[0] - legs[0]) * CUBE_RATE,
                    ((-1)**sequence[s][0] * o.c.size[0] + legs[0]) * CUBE_RATE],
                   [((-1)**sequence[s][1] * o.c.size[1] - legs[1]) * CUBE_RATE,
                    ((-1)**sequence[s][1] * o.c.size[1] + legs[1]) * CUBE_RATE],
                   [(-o.c.size[2] - legs[2]) * CUBE_RATE, (o.c.size[2] + legs[2]) * CUBE_RATE]] for s in range(n_legs)]
    # bound_cube = [[-o.c.size[0], o.c.size[1]], [-o.c.size[1], o.c.size[1]], [-o.c.size[2], o.c.size[2]]]

    for s in range(n_legs):
        # ind1 = 0 + int(s // 2 < 1)
        # ind2 = 1 + int(s // 2 < 2)
        r[0][(s + 1) * 6 + 0] = [bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0]]
        r[1][(s + 1) * 6 + 0] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 0] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 1] = [bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 1] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 1] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 2] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 2] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 2] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 3] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 3] = [bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 3] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 4] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 4] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 4] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0]]

        r[0][(s + 1) * 6 + 5] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 5] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 5] = [bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1]]

    for i in range(3):
        for k in range(2):
            tmp = k + i * 2
            ind1 = 0 + int(i < 1)
            ind2 = 1 + int(i < 2)
            for m in range(4):
                r[ind1][tmp] += [shift[ind1][sequence[m][0]]]
                r[ind2][tmp] += [shift[ind2][sequence[m][1]]]

    A = quart2dcm(o.c.q[j])
    for k in range(total_cubes):
        for i in range(4):
            r1 = A.T @ np.array([r[0][k][i], r[1][k][i], r[2][k][i]])
            r[0][k][i] = r1[0] + o.c.r_orf[j][0]
            r[1][k][i] = r1[1] + o.c.r_orf[j][1]
            r[2][k][i] = r1[2] + o.c.r_orf[j][2]
    anw = []
    for i in range(total_cubes):
        color = 'yellow' if i < 6 else 'gray'
        anw += [go.Mesh3d(x=r[0][i], y=r[1][i], z=r[2][i], color=color, opacity=1)]
    xl, yl, zl = ([], [], [])
    for i in range(int(len(o.c.line[j]) // 3)):
        xl += [o.c.line[j][i * 3 + 0]]
        yl += [o.c.line[j][i * 3 + 1]]
        zl += [o.c.line[j][i * 3 + 2]]
    anw += [go.Scatter3d(x=xl, y=yl, z=zl, mode='lines')]
    return anw

def show_chipsats_and_cubesats(o, clr: str = 'lightpink', opacity: float = 1):
    data = []
    for i in range(o.f.n):
        data += show_chipsat(o, i, clr, opacity)
    for i in range(o.c.n):
        data += show_cubesat(o, i)
    return go.Figure(data=data, layout=go.Layout(autosize=False, width=900, height=700))

def plot_all(o, save: bool = False, count: int = None):
    f = show_chipsats_and_cubesats(o, clr='black')
    # f.update_traces(scatter3d=dict(width=3))
    if save:
        f.write_image('img/' + str('{:04}'.format(count)) + '.jpg')
    else:
        f.show()

def plot_signals(o):
    global TITLE_SIZE, CAPTION_SIZE
    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            x = [o.p.show_rate * o.p.dt * i for i in range(len(o.c.signal_power[i_c][i_f]))]
            plt.plot(x, o.c.signal_power[i_c][i_f], label=f"фемтоступтник № {i_f+ 1}")
        plt.xlabel("Время, с", fontsize=CAPTION_SIZE)
        plt.ylabel("Мощность сигнала, ?", fontsize=CAPTION_SIZE)
        plt.title(f"График сигналов от фемтоспутников, получаемых кубсатом № {i_c + 1}")
        plt.legend(fontsize=CAPTION_SIZE)
        plt.show()
    tmp = plt.subplots(o.f.n, 1)
    fig = tmp[0]
    fig.suptitle(f"График сигналов от фемтоспутников, получаемых фемтосатом № {[i + 1 for i in range(o.f.n)]}",
                 fontsize=TITLE_SIZE)
    axes = tmp[1:o.f.n][0]
    colors = ['violet', 'teal', 'peru', 'cornflowerblue', 'forestgreen', 'blueviolet']
    for i_f1 in range(o.f.n):
        for i_f2 in range(o.f.n):
            if i_f1 != i_f2:
                x = [o.p.show_rate * o.p.dt * i for i in range(len(o.f.signal_power[i_f1][i_f2]))]
                axes[i_f1].plot(x, o.f.signal_power[i_f1][i_f2], c=colors[i_f2])
        axes[i_f1].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
        axes[i_f1].set_ylabel("?", fontsize=CAPTION_SIZE)
    plt.show()

def plot_distance(o):
    global TITLE_SIZE, CAPTION_SIZE
    fig, ax = plt.subplots(2, 2 if o.p.k.orientation else 1, figsize=(15 if o.p.k.orientation else 8, 10))
    axes = ax[0] if o.p.k.orientation else ax
    fig.suptitle(f"Неточности в навигации", fontsize=TITLE_SIZE)
    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            labels = ["Ошибка дистанции (реальная)",
                      "Ошибка дистанции (оцениваемая)",
                      "Ошибка определения положения"] \
                if i_f == 0 else [None for _ in range(100)]
            x = [o.p.dt * i for i in range(len(o.c.real_dist[i_c][i_f]))]
            axes[0].plot(x, np.abs(np.array(o.c.real_dist[i_c][i_f]) - np.array(o.c.calc_dist[i_c][i_f])),
                         c=MY_COLORS[0], label=labels[0] if i_c == 0 else None)
            axes[0].plot([o.p.dt * i for i in range(len(o.f.z_difference[i_f]))], o.f.z_difference[i_f],
                         c=MY_COLORS[3], label=labels[1] if i_c == 0 else None)
            axes[0].plot(x, [(-1)**i*10 if o.f.line_difference[i_f][i][0] == NO_LINE_FLAG else
                             np.linalg.norm(o.f.line_difference[i_f][i]) for i in range(len(x))],
                         c=MY_COLORS[2], label=labels[2] if i_c == 0 else None)
    axes[0].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
    axes[0].set_ylabel(f"Ошибка, м", fontsize=CAPTION_SIZE)
    axes[0].legend(fontsize=CAPTION_SIZE)
    axes[0].grid(True)

    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            labels = ["ΔX", "ΔY", "ΔZ"]
            x = [o.p.dt * i for i in range(len(o.c.real_dist[i_c][i_f]))]
            for j in range(3):
                axes[1].plot(x, [(-1)**i*10 if o.f.line_difference[i_f][i][j] == NO_LINE_FLAG else
                                 o.f.line_difference[i_f][i][j] for i in range(len(x))], c=MY_COLORS[j+3],
                             label=labels[j] if i_f == 0 and i_c == 0 else None)
    axes[1].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
    axes[1].set_ylabel(f"Компоненты r, м", fontsize=CAPTION_SIZE)
    axes[1].legend(fontsize=CAPTION_SIZE)
    axes[1].grid(True)

    if o.p.k.orientation:
        for i_f in range(o.f.n):
            labels_q = ["ΔΛ⁰", "ΔΛˣ", "ΔΛʸ", "ΔΛᶻ"]  # "ΔΛ⁰",
            labels_w = ["Δωˣ", "Δωʸ", "Δωᶻ"]  # "ΔΛ⁰",
            x = [o.p.dt * i for i in range(len(o.c.real_dist[0][i_f]))]
            for j in range(3):
                ax[1][0].plot(x, [o.f.attitude_difference[i_f][i][j] for i in range(len(x))],
                              c=MY_COLORS[j+3], label=labels_q[j] if i_f == 0 else None)
                ax[1][1].plot(x, [o.f.spin_difference[i_f][i][j] for i in range(len(x))],
                              c=MY_COLORS[j+3], label=labels_w[j] if i_f == 0 else None)
            ax[1][0].plot(x, [o.f.attitude_difference[i_f][i][3] for i in range(len(x))],
                          c=MY_COLORS[3+3], label=labels_q[3] if i_f == 0 else None)
        ax[1][0].set_ylabel(f"Компоненты Λ", fontsize=CAPTION_SIZE)
        ax[1][1].set_ylabel(f"Компоненты ω", fontsize=CAPTION_SIZE)
        for j in range(2):
            ax[1][j].legend(fontsize=CAPTION_SIZE)
            ax[1][j].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
            ax[1][j].grid(True)
    plt.show()

def plot_sigmas(o):
    """
    Функция берёт преобразованные диагональные элементы матрицы P и делает по ним выводы.
    Гордая такая.
    """
    fig, axes = plt.subplots(2, 1)
    fig.suptitle(f"Погрешности, оцениваемые фильтром Калмана", fontsize=TITLE_SIZE)
    colors = ['forestgreen', 'cornflowerblue', 'teal', 'peru']
    x = [o.p.dt * i for i in range(len(o.p.k.sigmas[0]))]
    t = 9 if o.p.k.orientation else 6
    for n in range(o.f.n):
        for i in range(t):
            j = int((i // 3) % 2 == 0)
            axes[j].plot(x, o.p.k.sigmas[n * t + i], c=colors[j])
            # axes[j].plot(x, o.p.k.real_sigmas[n * t + i], c=colors[j+2])
    '''for i in [0, 4]:
        j = int((i // 3) % 2 == 0)
        axes[j].plot(x, o.p.k.sigmas[i], c=colors[j], label="")
        # axes[j].plot(x, o.p.k.real_sigmas[i], c=colors[j+2])'''
    # plt.legend(fontsize=CAPTION_SIZE)
    plt.show()
