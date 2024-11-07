from PIL import Image
import plotly.graph_objs as go
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from my_math import *
from config import *


FEMTO_RATE = 1e2
CUBE_RATE = 1e2
TITLE_SIZE = 15  # 15
CAPTION_SIZE = 13  # 13
rcParams["savefig.directory"] = "/home/kodiak/Desktop"
rcParams["savefig.format"] = "jpg"

# >>>>>>>>>>>> 2D графики <<<<<<<<<<<<
def plot_distance(o):
    global TITLE_SIZE, CAPTION_SIZE
    fig, ax = plt.subplots(3 if o.v.NAVIGATION_ANGLES else 2,
                           2 if o.v.NAVIGATION_ANGLES else 1, figsize=(20 if o.v.NAVIGATION_ANGLES else 8, 10))
    axes = ax[0] if o.v.NAVIGATION_ANGLES else ax
    fig.suptitle(f"Неточности в навигации", fontsize=TITLE_SIZE)

    m = 0
    x = o.p.record['t'].to_list()

    for i_c in range(o.c.n):
        for i_f in range(o.f.n):
            labels = ["Ошибка дистанции (реальная)",
                      "Ошибка дистанции (оцениваемая)",
                      "Ошибка определения положения"] \
                if i_f == 0 else [None for _ in range(100)]
            y1 = o.p.record[f'{o.c.name}-{o.f.name} ErrorEstimateDistance {i_c} {i_f}'].to_list()
            y2 = o.p.record[f'ZModel&RealDifference'].to_list()
            y3 = o.p.record[f'{o.f.name} KalmanPosError r {i_f}'].to_list()
            axes[0].plot(x, y1, c=o.v.MY_COLORS[0], label=labels[0] if i_c == 0 else None)
            axes[0].plot(x, y2, c=o.v.MY_COLORS[3], label=labels[1] if i_c == 0 else None)
            axes[0].plot(x, y3, c=o.v.MY_COLORS[2], label=labels[2] if i_c == 0 else None)
            m = max(m, max(y1), max(y2), max(y3))
    axes[0].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
    axes[0].set_ylabel(f"Ошибка, м", fontsize=CAPTION_SIZE)
    axes[0].legend(fontsize=CAPTION_SIZE)
    axes[0].grid(True)
    if m > 1e3:
        axes[0].set_ylim([0, 1e3])

    for i_f in range(o.f.n):
        labels = ["ΔX", "ΔY", "ΔZ"]
        for j, c in enumerate('xyz'):
            y = o.p.record[f'{o.f.name} KalmanPosError {c} {i_f}'].to_list()
            axes[1].plot(x, y, c=o.v.MY_COLORS[j+3], label=labels[j] if i_f == 0 else None)
    axes[1].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
    axes[1].set_ylabel(f"Компоненты r, м", fontsize=CAPTION_SIZE)
    axes[1].legend(fontsize=CAPTION_SIZE)
    axes[1].grid(True)

    if o.v.NAVIGATION_ANGLES:
        for i_f in range(o.f.n):
            labels_q = ["Λˣ", "Λʸ", "Λᶻ"]
            labels_w = ["ωˣ", "ωʸ", "ωᶻ"]
            labels_dq = ["ΔΛˣ", "ΔΛʸ", "ΔΛᶻ"]
            labels_dw = ["Δωˣ", "Δωʸ", "Δωᶻ"]
            for j, c in enumerate('xyz'):
                y1 = o.p.record[f'{o.f.name} KalmanQuatError {c} {i_f}'].to_list()
                y2 = o.p.record[f'{o.f.name} KalmanSpinError ORF {c} {i_f}'].to_list()
                y3 = o.p.record[f'{o.f.name} RealQuat {c} {i_f}'].to_list()
                y3e = o.p.record[f'{o.f.name} KalmanQuatEstimation {c} {i_f}'].to_list()
                y4 = o.p.record[f'{o.f.name} RealSpin ORF {c} {i_f}'].to_list()
                y4e = o.p.record[f'{o.f.name} KalmanSpinEstimation ORF {c} {i_f}'].to_list()
                ax[1][0].plot(x, y1, c=o.v.MY_COLORS[j+3], label=labels_dq[j] if i_f == 0 else None)
                ax[1][1].plot(x, y2, c=o.v.MY_COLORS[j+3], label=labels_dw[j] if i_f == 0 else None)
                ax[2][0].plot(x, y3, label=labels_q[j] + "-реальное" if i_f == 0 else None)
                ax[2][1].plot(x, y4, label=labels_w[j] + "-реальное" if i_f == 0 else None)
                ax[2][0].plot(x, y3e, ":", label=labels_q[j] + "-оценка" if i_f == 0 else None)
                ax[2][1].plot(x, y4e, ":", label=labels_w[j] + "-оценка" if i_f == 0 else None)
        for ii in [1, 2]:
            ax[ii][0].set_ylabel(f"Компоненты Λ", fontsize=CAPTION_SIZE)
            ax[ii][1].set_ylabel(f"Компоненты ω (ORF)", fontsize=CAPTION_SIZE)
            for j in range(2):
                ax[ii][j].legend(fontsize=CAPTION_SIZE)
                ax[ii][j].set_xlabel("Время, с", fontsize=CAPTION_SIZE)
                ax[ii][j].grid(True)
    plt.show()

def plot_atmosphere_models(n: int = 100):
    from dynamics import get_atm_params
    v = Variables()

    range_km = [300, 500]
    fig, axs = plt.subplots(1, 2, figsize=(19, 5))
    fig.suptitle('Модели атмосферы')
    for m in range(len(v.ATMOSPHERE_MODELS)):
        for j in range(2):
            z = np.linspace(100e3 if j == 0 else range_km[0]*1e3, range_km[1]*1e3, n)
            rho = [get_atm_params(h=z[i], v=v, atm_model=v.ATMOSPHERE_MODELS[m])[0] for i in range(n)]
            tmp = ", используемая" if v.ATMOSPHERE_MODELS[m] == v.ATMOSPHERE_MODEL else ""
            axs[j].plot(z, rho, ls="-" if v.ATMOSPHERE_MODELS[m] == v.ATMOSPHERE_MODEL else ":",
                        label=f"Модель: {v.ATMOSPHERE_MODELS[m]}{tmp}")
            axs[j].set_ylabel(f"Плотность ρ, кг/м³")
            axs[j].set_xlabel(f"Высота z, м")
            axs[j].legend()
            axs[j].grid()
    axs[0].set_title(f"От линии Кармана до {range_km[1]} км")
    axs[1].set_title(f"От {range_km[0]} до {range_km[1]} км")
    plt.show()


# >>>>>>>>>>>> 3D отображение в ОСК <<<<<<<<<<<<
def show_chipsat(o, j, clr, opacity, reference_frame: str, return_go: bool = True, ax=None) -> list:
    """
    Функция отображения чипсата (просто пластина)
    :param o: Objects
    :param j: номер чипсата
    :param clr:
    :param opacity:
    :param reference_frame:
    :param return_go: По умолчанию, при return_go=False считается, что reference_frame="BRF"
    :param ax:
    :return:
    """
    global FEMTO_RATE
    rate = FEMTO_RATE if reference_frame != "BRF" else 2 / min(o.f.size)
    x, y, z = ([], [], [0 for _ in range(4)])
    for x_shift in [-o.f.size[0] * rate, o.f.size[0] * rate]:
        for y_shift in [-o.f.size[1] * rate, o.f.size[1] * rate]:
            x += [x_shift]
            y += [y_shift]
    A = quart2dcm(o.f.q[j])
    if reference_frame != "BRF":
        for i in range(4):
            r = A.T @ np.array([x[i], y[i], z[i]])  # НЕВЕРНО, что тут A.T !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            x[i] = r[0] + o.f.r_orf[j][0] if reference_frame == "ORF" else r[0] + o.f.r_irf[j][0]
            y[i] = r[1] + o.f.r_orf[j][1] if reference_frame == "ORF" else r[1] + o.f.r_irf[j][1]
            z[i] = r[2] + o.f.r_orf[j][2] if reference_frame == "ORF" else r[2] + o.f.r_irf[j][2]
        r_real, r_estimation = [], []
        for c in "xyz":
            r_real.append(o.p.record[f'{o.f.name} r {c} {reference_frame.lower()} {j}'].to_list())
            r_estimation.append(o.p.record[f'{o.f.name} KalmanPosEstimation {c} {j}'].to_list())
        if return_go:
            return [go.Mesh3d(x=x, y=y, z=z, color=clr, opacity=opacity),
                    go.Scatter3d(x=r_real[0], y=r_real[1], z=r_real[2], mode='lines'),
                    go.Scatter3d(x=r_estimation[0], y=r_estimation[1], z=r_estimation[2], mode='lines')]
    ax.plot([x[0], x[2], x[3], x[1], x[0]],
            [y[0], y[2], y[3], y[1], y[0]],
            [z[0], z[2], z[3], z[1], z[0]], c='gray', linewidth=3)

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
        anw.append(go.Mesh3d(x=r[0][i], y=r[1][i], z=r[2][i], color=color, opacity=1))
    r = []
    for c in 'xyz':
        r.append(o.p.record[f'{o.c.name} r {c} {reference_frame.lower()} {j}'].to_list())
    anw.append(go.Scatter3d(x=r[0], y=r[1], z=r[2], mode='lines'))
    return anw

def plot_model_gain(o: Objects, n: int = 20):
    from my_math import pol2dec
    from spacecrafts import get_gain
    fig = plt.figure(figsize=(15, 10))

    for i in range(2):
        for j in range(2):
            ax = fig.add_subplot(2, 2, i+1+j*2, projection='3d')
            obj = [o.f, o.c][i]
            my_print(f"Диаграмма направленностей для {obj.name}: {o.v.GAIN_MODEL_F}", color="b")

            u = np.linspace(0, 2 * np.pi, n)
            v = np.linspace(-np.pi / 2, np.pi / 2, n)
            U, V = np.meshgrid(u, v)

            max_g = 0
            for k in range(len(get_gain(v=o.v, obj=obj, r=pol2dec(1, u[0], v[0]), if_send=j == 0, if_take=j == 1))):
                g = np.array([[get_gain(v=o.v, obj=obj, r=pol2dec(1, u[ii], v[jj]),
                                           if_send=j == 0, if_take=j == 1)[k] for ii in range(n)] for jj in range(n)])
                X, Y, Z = pol2dec(g, U, V)
                ax.plot_surface(X, Y, Z, cmap=o.v.MY_COLORMAPS[k])
                max_g = max(max_g, np.max(g.flatten()))

            if obj.name == "FemtoSat":
                max_g = max(max_g, 2 * np.max(o.f.size))
                show_chipsat(o=o, j=0, reference_frame="BRF", return_go=False, ax=ax, clr=None, opacity=None)

            ax.set_xlim(-max_g, max_g)
            ax.set_ylim(-max_g, max_g)
            ax.set_zlim(-max_g, max_g)
            ax.set_box_aspect([1, 1, 1])
            title_text = f"Диаграмма направленностей для {obj.name} | " \
                         f"GAIN_MODEL = {o.v.GAIN_MODEL_F if  obj.name == 'FemtoSat' else o.v.GAIN_MODEL_C}\n" \
                         f"Отправленный сигнал"
            ax.set_title(title_text if j == 0 else "Принятый сигнал")
    plt.show()


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
    """Отрисовка слева красивой Земли из одной линии"""
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
    ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=bm, alpha=1)
    ax.set_xlabel("x, тыс. км")
    ax.set_ylabel("y, тыс. км")
    ax.set_zlabel("z, тыс. км")
    return ax

def plot_the_earth_go(v: Variables):
    spherical_earth_map = np.load('kiamformation/data/map_sphere.npy')
    xm, ym, zm = spherical_earth_map.T * v.EARTH_RADIUS

    return go.Scatter3d(x=xm, y=ym, z=zm, mode='lines')

def plot_reference_frames(ax, o, txt: str, color: str = "gray", t: float = None):
    from dynamics import get_matrices
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    z = np.array([0, 0, 1])
    arrows = np.array([x, y, z]) * o.v.ORBIT_RADIUS
    start = np.zeros(3)
    if txt == "ОСК":
        U, _, _, R_orb = get_matrices(obj=o.a, n=0, t=t, v=o.v)
        arrows = np.array([U.T @ x, U.T @ y, U.T @ z]) * o.v.ORBIT_RADIUS / 2
        start = R_orb
    if txt == "ИСК":
        # Отрисовка кружочка
        x, y, z = ([], [], [])
        """n_round = 50
        for t in np.linspace(0, o.v.SEC_IN_TURN * 1.5, n_round):
            _, _, _, R_orb = get_matrices(v=o.v, t=t, obj=o.a, n=0)"""
        for i in range(int(len(o.a.line_irf[0])//3)):
            R_orb = o.a.line_irf[0][3*i:3*i+3]
            x += [R_orb[0]]
            y += [R_orb[1]]
            z += [R_orb[2]]
        ax.plot(x, y, z, color)
    ax = arrows3d(starts=np.array([start for _ in range(3)]), ends=np.array([start + arrows[i] for i in range(3)]),
                  ax=ax, color=color, label=txt)
    for i in range(3):
        label = ["x", "y", "z"][i]
        a = start + arrows[i] + arrows[i] / np.linalg.norm(arrows[i]) * 0.2
        ax.text(a[0], a[1], a[2], c=color, s=label)
    return ax

# >>>>>>>>>>>> Анимация <<<<<<<<<<<<
def animate_reference_frames(resolution: int = 3, n: int = 5):
    from PIL import Image
    from os import remove

    v_ = Variables()
    o = Objects(v=v_)
    # o.v.dT = o.v.SEC_IN_TURN / (n - 3)
    TIME = 2*np.pi / o.v.W_ORB
    o.v.dT = TIME / n
    o.v.IF_NAVIGATION = False

    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlim3d([-1e7, 1e7])
    ax.set_ylim3d([-1e7, 1e7])
    ax.set_zlim3d([-1e7, 1e7])
    x_points = 180 * resolution
    y_points = 90 * resolution
    earth_image = Image.open(f'../../source/skins/{o.v.EARTH_FILE_NAME}')
    earth_image = np.array(earth_image.resize((x_points, y_points))) / 256.
    for i in range(n):
        o.p.time_step()
        ax = plot_the_earth_mpl(ax, v=v_, earth_image=earth_image)
        ax = plot_reference_frames(ax, o, t=o.p.t, txt="ИСК", color="lime")
        ax = plot_reference_frames(ax, o, t=o.p.t, txt="ОСК", color="red")
        ax.view_init(azim=20, elev=30, roll=0)
        ax.axis('equal')
        plt.title(f"Наклонение: {o.v.INCLINATION}°, эксцентриситет: {round(o.v.ECCENTRICITY, 2)}, "
                  f"апогей: {round(o.v.APOGEE / 1e3)} км, перигей: {round(o.v.PERIGEE / 1e3)} км")
        plt.legend()
        plt.savefig(f"../../res/to_delete_{'{:04}'.format(i)}.png")
        ax.clear()
    plt.close()

    images = [Image.open(f"../../res/to_delete_{'{:04}'.format(i)}.png") for i in range(n)]
    images[0].save('../../res/res.gif', save_all=True, append_images=images[1:], duration=20, loop=0)
    for i in range(n):
        remove(f"../../res/to_delete_{'{:04}'.format(i)}.png")

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
                                      f" | Время: {o.v.TIME} ({round(o.v.TIME / (3600 * 24), 2)} дней)  |  "
                                      f"i={o.v.INCLINATION}°, e={o.v.ECCENTRICITY}"))
    if save:
        fig.write_image('../../img/' + str('{:04}'.format(count)) + '.jpg')
    else:
        fig.show()
