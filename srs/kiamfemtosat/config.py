class Variables:
    def __init__(self):
        from kiam_astro import kiam
        import numpy
        from srs.kiamfemtosat.spacecrafts import Anchor
        # >>>>>>>>>>>> Вручную настраиваемые параметры <<<<<<<<<<<<
        self.dT = 1.
        self.TIME = 1e4
        self.CUBESAT_AMOUNT = 4
        self.CHIPSAT_AMOUNT = 1
        self.DYNAMIC_MODEL = {'aero drag': False,
                              'j2': False}
        self.NAVIGATION_BY_ALL = True
        self.NAVIGATION_ANGLES = False  # Содержит ли искомый вектор состояния кватернионы и угловые скорости
        self.MULTI_ANTENNA_TAKE = True  # Разделяет ли КА приходящий сигнал на составляющие
        self.MULTI_ANTENNA_SEND = True  # Разделяет ли КА исходящий сигнал на составляющие
        self.START_NAVIGATION_TOLERANCE = 0.9
        self.START_NAVIGATION = ['perfect', 'near', 'random'][1]
        self.GAIN_MODEL_C = ['isotropic', '1 antenna', '2 antennas', '3 antennas', 'ellipsoid'][0]
        self.GAIN_MODEL_F = ['isotropic', '1 antenna', '2 antennas', '3 antennas', 'ellipsoid'][0]
        self.SOLVER = ['rk4 hkw', 'kiamastro'][0]  # Везде проверяется на hkw -> проверки на rk4. Может изменить?
        self.CHIPSAT_OPERATING_MODE = ['const', 'while_sun_visible'][0]
        self.DISTORTION = 0.  # Искривление диаграммы направленности

        self.RVW_CubeSat_SPREAD = [1e2, 1e-1, 1e-5]  # r (м), v (м/с), ω (рад/с)
        self.RVW_ChipSat_SPREAD = [1e2, 1e-1, 1e-5]
        self.KALMAN_COEF = {'q': [1e-5]*2, 'p': [1e-5]*4, 'd': 1e-10, 'r': 1e-1}
        self.CUBESAT_MODEL = ['1U', '1.5U', '2U', '3U', '6U', '12U'][0]
        self.CHIPSAT_MODEL = ['KickSat', 'Трисат'][0]
        self.SHAMANISM = {'KalmanQuaternionNormalize': True,   # Нормировка кватернионов в фильтре Калмана
                          'KalmanSpinLimit': [True, 1e-2],  # Ограничение скорости вращения в прогнозе фильтра Калмана
                          'ClohessyWiltshireC1=0': True,  # Траектории без дрейфа (зануление C1, даже при аэродинамике)
                          'KalmanVelocityLimit': [True, 1e3],
                          'KalmanPositionLimit': [False, 1e4],
                          }

        self.ATMOSPHERE_MODEL = ['NASA', 'ПНБО', 'COESA62', 'COESA76'][0]  # Стояло 3 (20 сен)
        self.ATMOSPHERE_MODELS = ['NASA', 'ПНБО', 'COESA62', 'COESA76']

        self.N_ANTENNA_C = {'isotropic': 1, '1 antenna': 1, '2 antennas': 2, '3 antennas': 3,
                            'ellipsoid': 1}[self.GAIN_MODEL_C]
        self.N_ANTENNA_F = {'isotropic': 1, '1 antenna': 1, '2 antennas': 2, '3 antennas': 3,
                            'ellipsoid': 1}[self.GAIN_MODEL_F]


        # >>>>>>>>>>>> Параметры отображения <<<<<<<<<<<<
        self.IF_TALK = False
        self.IF_ANY_PRINT = True
        self.IF_TEST_PRINT = True
        self.IF_ANY_SHOW = False  # а ты к чему относишься?
        self.NO_LINE_FLAG = -10
        self.EARTH_FILE_NAME = ["earth1.jpg", "earth2.jpg", "earth3.webp"][2]


        # >>>>>>>>>>>> Константы <<<<<<<<<<<<
        self.ECCENTRICITY = 0.0
        self.INCLINATION = 0  # В градусах
        self.SEC_IN_TURN = 24*3600*kiam.units('earth')['TimeUnit']*2*numpy.pi
        self.SEC_IN_RAD = 24*3600*kiam.units('earth')['TimeUnit']
        self.EARTH_RADIUS = kiam.units('earth')['DistUnit'] * 1e3
        self.ORBIT_RADIUS = self.EARTH_RADIUS + 400e3

        # Параметры орбиты
        self.APOGEE = self.ORBIT_RADIUS  # Апогей
        self.PERIGEE = self.ORBIT_RADIUS * (1 - self.ECCENTRICITY)/(1 + self.ECCENTRICITY)  # Перигей
        self.P = self.APOGEE * (1 - self.ECCENTRICITY**2)  # Фокальный параметр
        self.MU = 5.972e24 * 6.67408e-11  # Гравитационный параметр
        self.W_ORB = numpy.sqrt(self.MU / self.ORBIT_RADIUS ** 3)
        self.W_ORB_VEC_IRF = self.W_ORB * numpy.array([0, -numpy.sin(self.INCLINATION), numpy.cos(self.INCLINATION)])
        self.V_ORB = numpy.sqrt(self.MU / self.ORBIT_RADIUS)
        self.J2 = 1.082 * 1e-3

        self.CUBESAT_MODELS = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        self.CHIPSAT_MODELS = ['KickSat', 'Трисат']
        self.GAIN_MODES = ['isotropic', '1 antenna', '2 antennas', '3 antennas', 'ellipsoid']
        self.NAVIGATIONS = ['perfect', 'near', 'random']
        self.OPERATING_MODES = ['free_flying', 'swarm_stabilize', 'lost']
        self.OPERATING_MODES_CHANGE = ['const', 'while_sun_visible']
        self.MY_COLORMAPS = ['cool', 'winter', 'summer', 'spring', 'gray', 'bone''autumn']
        self.MY_COLORS = ['violet', 'forestgreen', 'cornflowerblue', 'peru', 'teal', 'blueviolet', 'deeppink',
                          'darksalmon', 'magenta', 'maroon', 'orchid', 'purple', 'wheat', 'tan', 'steelblue',
                          'forestgreen',
                          'aqua', 'blue', 'beige', 'bisque', 'indigo', 'navy', 'deepskyblue', 'maroon', 'gold',
                          'aquamarine', 'indigo', 'olivedrab', 'slategray', 'pink', 'salmon', 'steelblue', 'peru']

        # >>>>>>>>>>>> Изменяемые параметры по ходу работы кода <<<<<<<<<<<<
        self.MEASURES_VECTOR = []
        self.MEASURES_VECTOR_NOTES = []

        # >>>>>>>>>>>> Параметры для тестов <<<<<<<<<<<<
        self.IF_NAVIGATION = True
        self.ANCHOR = Anchor(v=self)

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

    def plot_model_gain(n: int = 20):
        from srs.kiamfemtosat.my_plot import show_chipsat
        v_ = Variables()
        o = Objects(v=v_)
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
                for k in range(len(get_gain(v=v_, obj=obj, r=pol2dec(1, u[0], v[0]), if_send=j == 0, if_take=j == 1))):
                    g = np.array([[get_gain(v=v_, obj=obj, r=pol2dec(1, u[ii], v[jj]),
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

    def plot_atmosphere_models(n: int = 100):
        from srs.kiamfemtosat.dynamics import get_atm_params
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
                            label=f"Модель: {v.ATMOSPHERE_MODELS[m]}{tmp}")  # , c=MY_COLORS[m])
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

    # animate_reference_frames(resolution=1, n=30)
    plot_model_gain()
    # plot_atmosphere_models()
