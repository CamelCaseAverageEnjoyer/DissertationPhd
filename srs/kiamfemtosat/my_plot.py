from PIL import Image
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from matplotlib import animation
from srs.kiamfemtosat.my_math import *
from srs.kiamfemtosat.config import *


FEMTO_RATE = 1e2
CUBE_RATE = 1e2
TITLE_SIZE = 15  # 15
CAPTION_SIZE = 13  # 13
rcParams["savefig.directory"] = "/home/kodiak/Desktop"
rcParams["savefig.format"] = "jpg"

# >>>>>>>>>>>> 2D графики <<<<<<<<<<<<
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
                         c=o.v.MY_COLORS[0], label=labels[0] if i_c == 0 else None)
            axes[0].plot([o.p.dt * i for i in range(len(o.f.z_difference[i_f]))], o.f.z_difference[i_f],
                         c=o.v.MY_COLORS[3], label=labels[1] if i_c == 0 else None)
            axes[0].plot(x, [(-1)**i*10 if o.f.line_difference[i_f][i][0] == o.v.NO_LINE_FLAG else
                             np.linalg.norm(o.f.line_difference[i_f][i]) for i in range(len(x))],
                         c=o.v.MY_COLORS[2], label=labels[2] if i_c == 0 else None)
    axes[0].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
    axes[0].set_ylabel(f"Ошибка, м", fontsize=CAPTION_SIZE)
    axes[0].legend(fontsize=CAPTION_SIZE)
    axes[0].grid(True)

    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            labels = ["ΔX", "ΔY", "ΔZ"]
            x = [o.p.dt * i for i in range(len(o.c.real_dist[i_c][i_f]))]
            for j in range(3):
                axes[1].plot(x, [(-1)**i*10 if o.f.line_difference[i_f][i][j] == o.v.NO_LINE_FLAG else
                                 o.f.line_difference[i_f][i][j] for i in range(len(x))], c=o.v.MY_COLORS[j+3],
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
                              c=o.v.MY_COLORS[j+3], label=labels_q[j] if i_f == 0 else None)
                ax[1][1].plot(x, [o.f.spin_difference[i_f][i][j] for i in range(len(x))],
                              c=o.v.MY_COLORS[j+3], label=labels_w[j] if i_f == 0 else None)
            ax[1][0].plot(x, [o.f.attitude_difference[i_f][i][3] for i in range(len(x))],
                          c=o.v.MY_COLORS[3+3], label=labels_q[3] if i_f == 0 else None)
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


# >>>>>>>>>>>> 3D отображение в ОСК <<<<<<<<<<<<
def show_chipsat(o, j, clr, opacity, reference_frame: str) -> list:
    global FEMTO_RATE
    x, y, z = ([], [], [0 for _ in range(4)])
    for x_shift in [-o.f.size[0] * FEMTO_RATE, o.f.size[0] * FEMTO_RATE]:
        for y_shift in [-o.f.size[1] * FEMTO_RATE, o.f.size[1] * FEMTO_RATE]:
            x += [x_shift]
            y += [y_shift]
    A = quart2dcm(o.f.q[j])
    for i in range(4):
        r = A.T @ np.array([x[i], y[i], z[i]])
        x[i] = r[0] + o.f.r_orf[j][0] if reference_frame == "ORF" else r[0] + o.f.r_irf[j][0]
        y[i] = r[1] + o.f.r_orf[j][1] if reference_frame == "ORF" else r[1] + o.f.r_irf[j][1]
        z[i] = r[2] + o.f.r_orf[j][2] if reference_frame == "ORF" else r[2] + o.f.r_irf[j][2]
    xl, yl, zl = ([], [], [])
    for i in range(int(len(o.f.line_orf[j]) // 3)):
        xl += [o.f.line_orf[j][i * 3 + 0]] if reference_frame == "ORF" else [o.f.line_irf[j][i * 3 + 0]]
        yl += [o.f.line_orf[j][i * 3 + 1]] if reference_frame == "ORF" else [o.f.line_irf[j][i * 3 + 1]]
        zl += [o.f.line_orf[j][i * 3 + 2]] if reference_frame == "ORF" else [o.f.line_irf[j][i * 3 + 2]]
    xk, yk, zk = ([], [], [])
    for i in range(int(len(o.f.line_kalman[j])//3)):
        xk += [o.f.line_kalman[j][i * 3 + 0]]
        yk += [o.f.line_kalman[j][i * 3 + 1]]
        zk += [o.f.line_kalman[j][i * 3 + 2]]
    return [go.Mesh3d(x=x, y=y, z=z, color=clr, opacity=opacity), go.Scatter3d(x=xl, y=yl, z=zl, mode='lines'),
            go.Scatter3d(x=xk, y=yk, z=zk, mode='lines')]

def show_cubesat(o, j, reference_frame: str) -> list:
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
            r1 = A.T @ np.array([r[0][k][i], r[1][k][i], r[2][k][i]])  # При IRF нет поворота. А и хрен с ним
            r[0][k][i] = r1[0] + o.c.r_orf[j][0] if reference_frame == "ORF" else r1[0] + o.c.r_irf[j][0]
            r[1][k][i] = r1[1] + o.c.r_orf[j][1] if reference_frame == "ORF" else r1[1] + o.c.r_irf[j][1]
            r[2][k][i] = r1[2] + o.c.r_orf[j][2] if reference_frame == "ORF" else r1[2] + o.c.r_irf[j][2]
    anw = []
    for i in range(total_cubes):
        color = 'yellow' if i < 6 else 'gray'
        anw += [go.Mesh3d(x=r[0][i], y=r[1][i], z=r[2][i], color=color, opacity=1)]
    xl, yl, zl = ([], [], [])
    for i in range(int(len(o.c.line_orf[j]) // 3)):
        xl += [o.c.line_orf[j][i * 3 + 0]] if reference_frame == "ORF" else [o.c.line_irf[j][i * 3 + 0]]
        yl += [o.c.line_orf[j][i * 3 + 1]] if reference_frame == "ORF" else [o.c.line_irf[j][i * 3 + 1]]
        zl += [o.c.line_orf[j][i * 3 + 2]] if reference_frame == "ORF" else [o.c.line_irf[j][i * 3 + 2]]
    anw += [go.Scatter3d(x=xl, y=yl, z=zl, mode='lines')]
    return anw


# >>>>>>>>>>>> 3D отображение в ИСК <<<<<<<<<<<<
def arrows3d(ends: np.ndarray, starts: np.ndarray = None, ax=None, label: str = None, **kwargs):
    """Построение 3D стрелок
    GitHub: https://github.com/matplotlib/matplotlib/issues/22571
    :param ends: (N, 3) size array of arrow end coordinates
    :param starts: (N, 3) size array of arrow start coordinates.
    :param ax: (Axes3DSubplot) existing axes to add to
    :param label: legend label to apply to this group of arrows
    :param kwargs: additional arrow properties"""
    if starts is None:
        starts = np.zeros_like(ends)

    assert starts.shape == ends.shape, "`starts` and `ends` shape must match"
    assert len(ends.shape) == 2 and ends.shape[1] == 3, \
        "`starts` and `ends` must be shape (N, 3)"

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            super().__init__((0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def do_3d_projection(self, renderer=None):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

            return np.min(zs)

    # create new axes if none given
    if ax is None:
        ax = plt.figure().add_subplot(111, projection='3d')
    arrow_prop_dict = dict(mutation_scale=20, arrowstyle='-|>', color='k', shrinkA=0, shrinkB=0)
    arrow_prop_dict.update(kwargs)
    for ind, (s, e) in enumerate(np.stack((starts, ends), axis=1)):
        a = Arrow3D(
            [s[0], e[0]], [s[1], e[1]], [s[2], e[2]],
            # only give label to first arrow
            label=label if ind == 0 else None,
            **arrow_prop_dict)
        ax.add_artist(a)
    ax.points = np.vstack((starts, ends, getattr(ax, 'points', np.empty((0, 3)))))
    return ax

def plot_the_earth_mpl(ax, v: Variables, res: int = 1, pth: str = "./", earth_image=None):
    x_points = 180 * res
    y_points = 90 * res

    if earth_image is None:
        bm = Image.open(f'{pth}source/skins/{v.EARTH_FILE_NAME}')
        bm = np.array(bm.resize((x_points, y_points))) / 256.
    else:
        bm = earth_image

    lons = np.linspace(-180, 180, bm.shape[1]) * np.pi / 180
    lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi / 180

    x = v.EARTH_RADIUS * np.outer(np.cos(lons), np.cos(lats)).T
    y = v.EARTH_RADIUS * np.outer(np.sin(lons), np.cos(lats)).T
    z = v.EARTH_RADIUS * np.outer(np.ones(np.size(lons)), np.sin(lats)).T
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=bm)
    ax.set_xlabel("x, км")
    ax.set_ylabel("y, км")
    ax.set_zlabel("z, км")
    return ax

def plot_the_earth_go(v: Variables):
    spherical_earth_map = np.load('data/map_sphere.npy')
    xm, ym, zm = spherical_earth_map.T * v.EARTH_RADIUS

    return go.Scatter3d(x=xm, y=ym, z=zm, mode='lines')

def plot_reference_frames(ax, o, txt: str, color: str = "gray", t: float = None):
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    arrows = np.array([x, y, z]) * o.v.ORBIT_RADIUS
    start = np.zeros(3)
    if txt == "ОСК":
        # Костыль на отсутствие вращения тела
        # get_matrices(self, obj: Union[CubeSat, FemtoSat], n: int, t: float = None)
        o.c.q[0] = [1, 0, 0, 0]
        U, S, A, R_orb = o.p.get_matrices(obj=o.c, n=0, t=t)
        arrows = np.array([U.T @ x, U.T @ y, U.T @ z]) * o.v.ORBIT_RADIUS / 2
        start = R_orb
    if txt == "ИСК":
        n_round = 30
        ax.plot(o.v.ORBIT_RADIUS * np.array([np.cos(i) for i in np.linspace(0, 2 * np.pi, n_round)]),
                o.v.ORBIT_RADIUS * np.array([np.sin(i) for i in np.linspace(0, 2 * np.pi, n_round)]), np.zeros(n_round),
                color)  # , ls=":"
    ax = arrows3d(starts=np.array([start for _ in range(3)]), ends=np.array([start + arrows[i] for i in range(3)]),
                  ax=ax, color=color, label=txt)
    for i in range(3):
        label = ["x", "y", "z"][i]
        a = start + arrows[i] + arrows[i] / np.linalg.norm(arrows[i]) * 0.2
        ax.text(a[0], a[1], a[2], c=color, s=label)
    return ax


# >>>>>>>>>>>> Конечная ассамблея <<<<<<<<<<<<
def show_chipsats_and_cubesats(o, reference_frame: str, clr: str = 'lightpink', opacity: float = 1):
    data = []
    for i in range(o.f.n):
        data += show_chipsat(o, i, clr, opacity, reference_frame)
    for i in range(o.c.n):
        data += show_cubesat(o, i, reference_frame)
    return data

def plot_all(o, save: bool = False, count: int = None):
    from plotly.subplots import make_subplots
    fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'surface'}, {'type': 'surface'}]],
                        subplot_titles=('Инерциальная СК', 'Орбитальная СК'))
    for i in range(2):
        tmp = show_chipsats_and_cubesats(o, ['IRF', 'ORF'][i], clr='black')
        for surf in tmp:
            fig.add_trace(surf, row=1, col=i+1)
    fig.add_trace(plot_the_earth_go(v=o.v), row=1, col=1)
    fig.update_layout(title=dict(text=f"Солвер: {o.v.SOLVER} "
                                      f"{'(' if o.v.DYNAMIC_MODEL['aero drag'] or o.v.DYNAMIC_MODEL['j2'] else ''}"
                                      f"{' +Лобовое сопротивление' if o.v.DYNAMIC_MODEL['aero drag'] else ''}"
                                      f"{' +Вторая гармоника' if o.v.DYNAMIC_MODEL['j2'] else ''}"
                                      f"{' )' if o.v.DYNAMIC_MODEL['aero drag'] or o.v.DYNAMIC_MODEL['j2'] else ''}"
                                      f" | Время: {o.v.TIME} (дней: {round(o.v.TIME / (3600 * 24), 2)})"))
    if save:
        fig.write_image('img/' + str('{:04}'.format(count)) + '.jpg')
    else:
        fig.show()
