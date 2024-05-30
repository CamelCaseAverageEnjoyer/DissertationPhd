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
GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas']
NAVIGATIONS = ['perfect', 'near', 'random']
OPERATING_MODES = ['free_flying', 'swarm_stabilize', 'lost']
OPERATING_MODES_CHANGE = ['const', 'while_sun_visible']
MY_COLORMAPS = ['cool', 'winter', 'summer', 'spring', 'gray', 'bone''autumn', ]
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

    plot_model_gain()
