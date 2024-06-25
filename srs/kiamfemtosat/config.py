class Variables:
    def __init__(self):
        import numpy
        # >>>>>>>>>>>> Вручную настраиваемые параметры <<<<<<<<<<<<
        self.dT = 10.
        self.TIME = 1e4
        self.CUBESAT_AMOUNT = 1
        self.CHIPSAT_AMOUNT = 1
        self.AERO_DRAG = False
        self.NAVIGATION_BY_ALL = True
        self.NAVIGATION_ANGLES = False
        self.START_NAVIGATION_TOLERANCE = 0.9
        self.DYNAMIC_MODEL = ['rk4', 'kiamastro'][1]
        self.START_NAVIGATION = ['perfect', 'near', 'random'][2]
        self.GAIN_MODEL_C = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas'][5]
        self.GAIN_MODEL_F = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas'][4]
        self.CHIPSAT_OPERATING_MODE = ['const', 'while_sun_visible'][0]
        self.DISTORTION = 0.  # Искривление диаграммы направленности

        self.RVW_CubeSat_SPREAD = [1e1, 1e-1, 1e-5]
        self.RVW_ChipSat_SPREAD = [1e2, 1e-1, 1e-5]
        self.KALMAN_COEF = {'q': [1e-12]*2, 'p': [1e-7]*4, 'r': 1e-1}
        self.CUBESAT_MODEL = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
        self.SHAMANISM = {'KalmanQuaternionNormalize': True,   # Нормировка кватернионов в фильтре Калмана
                          'KalmanSpinLimit': [True, 1e-3],  # Ограничение скорости вращения в прогнозе фильтра Калмана
                          'ClohessyWiltshireC1=0': True}  # Траектории без дрейфа (зануление C1, даже при аэродинамике)

        self.ATMOSPHERE_MODEL = ['NASA', 'ПНБО', 'COESA62', 'COESA76'][3]


        # >>>>>>>>>>>> Параметры отображения <<<<<<<<<<<<
        self.IF_TALK = False
        self.IF_ANY_PRINT = True
        self.NO_LINE_FLAG = -10
        self.EARTH_FILE_NAME = ["earth1.jpg", "earth2.jpg", "earth3.webp"][2]


        # >>>>>>>>>>>> Константы <<<<<<<<<<<<
        self.EARTH_RADIUS = 6731e3
        self.ORBIT_RADIUS = self.EARTH_RADIUS + 400e3
        self.ECCENTRICITY = 0
        self.INCLINATION = 0  # В градусах

        self.MU = 5.972e24 * 6.67408e-11  # гравитационный параметр
        self.W_ORB = numpy.sqrt(self.MU / self.ORBIT_RADIUS ** 3)
        self.V_ORB = numpy.sqrt(self.MU / self.ORBIT_RADIUS)
        self.J2 = 1.082 * 1e-3
        self.CUBESAT_MODELS = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        self.GAIN_MODES = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas', '1+1+1 antennas']
        self.NAVIGATIONS = ['perfect', 'near', 'random']
        self.OPERATING_MODES = ['free_flying', 'swarm_stabilize', 'lost']
        self.OPERATING_MODES_CHANGE = ['const', 'while_sun_visible']
        self.MY_COLORMAPS = ['cool', 'winter', 'summer', 'spring', 'gray', 'bone''autumn']
        self.MY_COLORS = ['violet', 'forestgreen', 'cornflowerblue', 'peru', 'teal', 'blueviolet', 'deeppink',
                          'darksalmon', 'magenta', 'maroon', 'orchid', 'purple', 'wheat', 'tan', 'steelblue',
                          'forestgreen',
                          'aqua', 'blue', 'beige', 'bisque', 'indigo', 'navy', 'deepskyblue', 'maroon', 'gold',
                          'aquamarine', 'indigo', 'olivedrab', 'slategray', 'pink', 'salmon', 'steelblue', 'peru']


        # >>>>>>>>>>>> Параметры для тестов <<<<<<<<<<<<
        self.IF_NAVIGATION = False

    def test_mode(self):
        self.IF_TALK = False
        self.IF_ANY_PRINT = False
        self.IF_NAVIGATION = False


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from srs.kiamfemtosat.spacecrafts import *
    from srs.kiamfemtosat.my_math import pol2dec
    from srs.kiamfemtosat.main import Objects
    from srs.kiamfemtosat.cosmetic import my_print

    def plot_model_gain(n: int = 30):
        v_ = Variables()
        o = Objects(v=v_)
        my_print(f"Диаграмма направленностей для ChipSat: {o.v.GAIN_MODEL_F}", color="b")
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection='3d')

        u = np.linspace(0, 2 * np.pi, n)
        v = np.linspace(-np.pi / 2, np.pi / 2, n)
        U, V = np.meshgrid(u, v)

        max_g = 0
        for k in range(len(get_gain(v=v_, obj=o.f, r=pol2dec(1, u[0], v[0])))):
            g = np.array([[get_gain(v=v_, obj=o.f, r=pol2dec(1, u[i], v[j]))[k] for i in range(n)] for j in range(n)])
            X, Y, Z = pol2dec(g, U, V)
            ax.plot_surface(X, Y, Z, cmap=o.v.MY_COLORMAPS[k])
            max_g = max(max_g, np.max(g.flatten()))

        ax.set_xlim(-max_g, max_g)
        ax.set_ylim(-max_g, max_g)
        ax.set_zlim(-max_g, max_g)
        ax.set_box_aspect([1, 1, 1])
        ax.set_title(f"Диаграмма направленностей для ChipSat | GAIN_MODEL = {o.v.GAIN_MODEL_F}")
        plt.show()

    def plot_atmosphere_models(n: int = 100):
        from srs.kiamfemtosat.dynamics import get_atm_params
        atmosphere_models = ['NASA', 'ПНБО', 'COESA62', 'COESA76']
        v = Variables()

        range_km = [300, 500]
        fig, axs = plt.subplots(1, 2, figsize=(19, 5))
        fig.suptitle('Модели атмосферы')
        for m in range(len(atmosphere_models)):
            for j in range(2):
                z = np.linspace(100e3 if j == 0 else range_km[0]*1e3, range_km[1]*1e3, n)
                rho = [get_atm_params(h=z[i], v=v, atm_model=atmosphere_models[m])[0] for i in range(n)]
                tmp = ", используемая" if atmosphere_models[m] == v.ATMOSPHERE_MODEL else ""
                axs[j].plot(z, rho, ls="-" if atmosphere_models[m] == v.ATMOSPHERE_MODEL else ":",
                            label=f"Модель: {atmosphere_models[m]}{tmp}")  # , c=MY_COLORS[m])
                axs[j].set_ylabel(f"Плотность ρ, кг/м³")
                axs[j].set_xlabel(f"Высота z, м")
                axs[j].legend()
                axs[j].grid()
        axs[0].set_title(f"От линии Кармана до {range_km[1]} км")
        axs[1].set_title(f"От {range_km[0]} до {range_km[1]} км")
        plt.show()

    def animate_reference_frames(resolution: int = 3, n: int = 5):
        from srs.kiamfemtosat.my_plot import plot_the_earth_mpl, plot_reference_frames
        from srs.kiamfemtosat.main import Objects
        from PIL import Image
        from os import remove

        v_ = Variables()
        o = Objects(v=v_)
        o.v.dT = 1.5 * 3600 / n
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
            # t = t_list[i] * 3600  # t в секундах, i в часах
            ax = plot_the_earth_mpl(ax, v=v_, earth_image=earth_image)
            ax = plot_reference_frames(ax, o, txt="ИСК", color="lime")
            # ax = plot_reference_frames(ax, o, txt="ОСК", color="red", t=t)
            ax = plot_reference_frames(ax, o, txt="ОСК", color="red")
            ax.axis('equal')
            plt.legend()
            plt.savefig(f"../../res/to_delete_{'{:04}'.format(i)}.jpg")
            ax.clear()
        plt.close()

        images = [Image.open(f"../../res/to_delete_{'{:04}'.format(i)}.jpg") for i in range(n)]
        images[0].save('../../res/res.gif', save_all=True, append_images=images[1:], duration=100, loop=0)
        for i in range(n):
            remove(f"../../res/to_delete_{'{:04}'.format(i)}.jpg")

    # animate_reference_frames(resolution=1, n=15)
    plot_model_gain()
    plot_atmosphere_models()
