from srs.kiamfemtosat.dynamics import *

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
                             kalman_coef=KALMAN_COEF if kalman_coef is None else kalman_coef, h_orb=self.h_orb)

    def integrate(self, t: float, animate: bool = False) -> None:
        my_print(f"Оборотов вокруг Земли: {round(10 * t / (3600 * 1.5)) / 10}  "
                 f"(дней: {round(100 * t / (3600 * 24)) / 100})", color='b')
        n = int(t // self.p.dt)
        flag = [0., 0.]
        for i in range(n):
            if i / n > (flag[0] + 0.1):
                flag[0] += 0.1
                per = int(10 * i / n)
                my_print(f"{10 * per}% [{'#' * per + ' ' * (10 - per)}]" +
                         real_workload_time(n=per, n_total=10, time_begin=self.p.time_begin,
                                            time_now=datetime.now()), color='m')
            if animate and i / n > (flag[1] + 0.01):
                flag[1] += 0.01
                plot_all(self, save=True, count=int(flag[1] // 0.01))
            self.p.time_step()

if __name__ == '__main__':
    # Инициализация
    o = Objects(n_c=CUBESAT_AMOUNT, n_f=CHIPSAT_AMOUNT, dt=dT, start_navigation_tolerance=0.9,
                navigation_method='kalman_filter rv' + 'q' * int(NAVIGATION_ANGLES) + ' all' * int(NAVIGATION_BY_ALL),
                kalman_coef=KALMAN_COEF, model_c=CUBESAT_MODEL, start_navigation=START_NAVIGATION)
    o.f.gain_mode = GAIN_MODEL
    o.c.gain_mode = GAIN_MODEL
    o.p.k.gain_mode = o.f.gain_mode
    o.p.is_aero = AERO_DRAG
    o.p.show_rate = 1

    # Интегрирование
    o.integrate(t=TIME, animate=False)

    # Вывод результатов
    tmp = np.array([np.linalg.norm(o.f.line_difference[0][i]) for i in range(len(o.f.line_difference[0]))])
    print(f"Математическое ожидание: {tmp.mean()}, Среднее отклонение: {tmp.std()}")
    # talk_decision()
    # plot_all(o)
    plot_distance(o)
    # plot_sigmas(o)
