import numpy as np

# >>>>>>>>>>>> Вручную настраиваемые параметры <<<<<<<<<<<<
dT = 10.
TIME = 1e5
CUBESAT_AMOUNT = 5
CHIPSAT_AMOUNT = 2
AERO_DRAG = False
NAVIGATION_BY_ALL = True
NAVIGATION_ANGLES = False
START_NAVIGATION = ['perfect', 'near', 'random'][2]
GAIN_MODEL = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas'][0]
DISTORTION = 0.  # Искривление диаграммы направленности
CHIPSAT_OPERATING_MODE = ['const', 'while_sun_visible'][0]

SHAMANISM = {'KalmanQuaternionNormalize': True,   # Нормировка кватернионов в фильтре Калмана
             'KalmanSpinLimit': [True, 1e-3],  # Ограничение скорости вращения в прогнозе фильтра Калмана
             'ClohessyWiltshireC1=0': True}  # Траектории без дрейфа (зануление C1, при учёте аэродинамики поломок нет)
KALMAN_COEF = {'q': [1e-12]*2, 'p': [1e-7]*4, 'r': 1e-1}
R_V_CubeSat_SPREAD = [10, 0.1]
CUBESAT_MODEL = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]

# >>>>>>>>>>>> Параметры отображения <<<<<<<<<<<<
IF_ANY_PRINT = True
NO_LINE_FLAG = -10

# >>>>>>>>>>>> Константы <<<<<<<<<<<<
EARTH_RADIUS = 6731e3
ORBIT_RADIUS = EARTH_RADIUS + 400e3
MU = 5.972e24 * 6.67408e-11  # гравитационный параметр
W_ORB = np.sqrt(MU / ORBIT_RADIUS ** 3)
V_ORB = np.sqrt(MU / ORBIT_RADIUS)
J2 = 1.082 * 1e-3
EARTH_FILE_NAME = ["earth1.jpg", "earth2.jpg", "earth3.webp"][2]
GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas']
NAVIGATIONS = ['perfect', 'near', 'random']
OPERATING_MODES = ['free_flying', 'swarm_stabilize', 'lost']
OPERATING_MODES_CHANGE = ['const', 'while_sun_visible']
MY_COLORMAPS = ['cool', 'winter', 'summer', 'spring', 'gray', 'bone''autumn']
MY_COLORS = ['violet', 'blueviolet', 'forestgreen', 'cornflowerblue', 'peru', 'teal', 'blueviolet', 'deeppink',
             'darksalmon', 'magenta', 'maroon', 'orchid', 'purple', 'wheat', 'tan', 'steelblue', 'forestgreen',
             'aqua', 'blue', 'beige', 'bisque', 'indigo', 'navy', 'deepskyblue', 'maroon', 'gold', 'aquamarine',
             'indigo', 'olivedrab', 'slategray', 'pink', 'salmon', 'steelblue', 'peru']


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    from srs.kiamfemtosat.spacecrafts import *
    from srs.kiamfemtosat.my_math import pol2dec
    from srs.kiamfemtosat.main import Objects
    from srs.kiamfemtosat.cosmetic import my_print


    def plot_model_gain(N: int = 30):
        my_print(f"Диаграмма направленностей типа: {GAIN_MODEL}", color="b")
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection='3d')

        u = np.linspace(0, 2 * np.pi, N)
        v = np.linspace(-np.pi / 2, np.pi / 2, N)
        U, V = np.meshgrid(u, v)

        o_ = Objects()
        o_.f.gain_mode = GAIN_MODEL
        max_g = 0

        for k in range(len(get_gain(obj=o_.f, r=pol2dec(1, u[0], v[0])))):
            g = np.array([[get_gain(obj=o_.f, r=pol2dec(1, u[i], v[j]))[k] for i in range(N)] for j in range(N)])
            X, Y, Z = pol2dec(g, U, V)
            ax.plot_surface(X, Y, Z, cmap=MY_COLORMAPS[k])
            max_g = max(max_g, np.max(g.flatten()))

        ax.set_xlim(-max_g, max_g)
        ax.set_ylim(-max_g, max_g)
        ax.set_zlim(-max_g, max_g)
        ax.set_box_aspect([1, 1, 1])
        ax.set_title(f"Диаграмма направленностей | GAIN_MODEL = {GAIN_MODEL}")
        plt.show()

    def animate_reference_frames(resolution: int = 3, n: int = 5):
        from srs.kiamfemtosat.my_plot import plot_the_earth, plot_reference_frames
        from srs.kiamfemtosat.main import Objects
        from PIL import Image
        from os import remove

        o = Objects()
        t_list = np.linspace(0, 1.5, n)
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(projection='3d')
        ax.set_xlim3d([-1e7, 1e7])
        ax.set_ylim3d([-1e7, 1e7])
        ax.set_zlim3d([-1e7, 1e7])
        x_points = 180 * resolution
        y_points = 90 * resolution
        earth_image = Image.open(f'../../source/skins/{EARTH_FILE_NAME}')
        earth_image = np.array(earth_image.resize((x_points, y_points))) / 256.
        for i in range(n):
            t = t_list[i] * 3600  # t в секундах, i в часах
            ax = plot_the_earth(ax, earth_image=earth_image)
            ax = plot_reference_frames(ax, o, txt="ИСК", color="lime")
            ax = plot_reference_frames(ax, o, txt="ОСК", color="red", t=t)
            ax.axis('equal')
            plt.legend()
            plt.savefig(f"../../res/to_delete_{'{:04}'.format(i)}.jpg")
            ax.clear()
        plt.close()

        images = [Image.open(f"../../res/to_delete_{'{:04}'.format(i)}.jpg") for i in range(n)]
        images[0].save('../../res/res.gif', save_all=True, append_images=images[1:], duration=100, loop=0)
        for i in range(n):
            remove(f"../../res/to_delete_{'{:04}'.format(i)}.jpg")

    # animate_reference_frames(resolution=5, n=15)
    plot_model_gain()
