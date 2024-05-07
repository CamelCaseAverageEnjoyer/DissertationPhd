from dynamics import *
# from navigation import *
# from guidance import *
# from control import *


class Objects:
    def __init__(self, navigation_method: str = None, kalman_coef: dict = None, dt: float = 1,
                 n_c: int = 1, n_f: int = 5, model_c: str = '1U',
                 start_navigation: str = 'perfect', start_navigation_tolerance: float = 0.9):
        """Класс содержит информацию о n_c кубсатах и n_f фемтосатах. \n
        Размер n_c кубсатов определяется моделью model_c"""
        # Проверки на вшивость
        models = ['1U', '1.5U', '2U', '3U', '6U', '12U']
        if model_c not in models:
            raise ValueError("Внимание! Модель кубсата должна быть одна из: 1U/1.5U/2U/3U/6U/12U")

        self.n_a = 1
        self.n_f = n_f
        self.model_c = model_c

        # Классы
        self.h_orb = 600e3
        self.mu = 5.972e24 * 6.67408e-11  # гравитационный параметр
        self.j_2 = 1.082 * 1e-3
        self.r_earth = 6371e3
        self.r_orb = self.r_earth + self.h_orb
        self.w_orb = np.sqrt(self.mu / self.r_orb ** 3)
        self.c = CubeSat(n=n_c, n_f=n_f, model=model_c, w_orb=self.w_orb)
        self.f = FemtoSat(n=n_f, start_navigation_tolerance=start_navigation_tolerance,
                          start_navigation=start_navigation, w_orb=self.w_orb)
        self.p = PhysicModel(c=self.c, f=self.f, dt=dt,
                             method_navigation='kalman_filter rv' if navigation_method is None else navigation_method,
                             kalman_coef={'q': 1e-12, 'p': [1e-8] * 3, 'r': 1e-1}
                             if kalman_coef is None else kalman_coef, h_orb=self.h_orb)


if __name__ == '__main__':
    # Инициализация
    t_integrate = 1e5
    o = Objects(n_c=1, n_f=1, dt=1., start_navigation_tolerance=0.9,
                navigation_method='kalman_filter rv' + 'q' * 0 + ' all' * 0,
                kalman_coef={'q': [1e-10, 1e-15], 'p': [1e-8]*4, 'r': 1e1},
                model_c=['1U', '1.5U', '2U', '3U', '6U', '12U'][0],
                start_navigation=['perfect', 'near', 'random'][1])
    o.f.gain_mode = ['isotropic', 'ellipsoid', '1 antenna', '2 antennas', '1+1 antennas'][0]
    o.c.gain_mode = o.f.gain_mode
    o.p.k.gain_mode = o.f.gain_mode
    o.p.is_aero = False
    o.p.show_rate = 1

    # Интегрирование
    o.p.integrate(t=t_integrate, animate=False)

    # Результаты
    tmp = np.array([np.linalg.norm(o.f.line_difference[0][i]) for i in range(len(o.f.line_difference[0]))])
    print(f"Математическое ожидание: {tmp.mean()}, Среднее отклонение: {tmp.std()}")
    # plot_all(o)
    plot_distance(o)
    # plot_sigmas(o)
