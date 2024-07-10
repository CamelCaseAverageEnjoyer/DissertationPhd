"""Численное моделирование космической миссии с использованием чипсатов"""
from srs.kiamfemtosat.dynamics import *
from srs.kiamfemtosat.cosmetic import *

class Objects:
    def __init__(self, v: Variables):
        """Класс объединяет следующие другие классы: CubeSat, FemtoSat, PhysicModel"""
        # Классы
        self.v = v
        self.a = v.ANCHOR
        self.c = CubeSat(v=v)
        self.f = FemtoSat(v=v)
        self.p = PhysicModel(c=self.c, f=self.f, a=self.a, v=v)

    def integrate(self, t: float, animate: bool = False) -> None:
        my_print(f"Оборотов вокруг Земли: {round(t / (2 * np.pi / self.v.W_ORB), 2)}  "
                 f"(дней: {round(t / (3600 * 24), 2)})", color='b', if_print=self.v.IF_ANY_PRINT)
        n = int(t // self.v.dT)
        flag = [0., 0.]
        for i in range(n):
            if i / n > (flag[0] + 0.1):
                flag[0] += 0.1
                per = int(10 * i / n)
                my_print(f"{10 * per}% [{'#' * per + ' ' * (10 - per)}]" +
                         real_workload_time(n=per, n_total=10, time_begin=self.p.time_begin,
                                            time_now=datetime.now()), color='m', if_print=self.v.IF_ANY_PRINT)
            if animate and i / n > (flag[1] + 0.01):
                flag[1] += 0.01
                plot_all(self, save=True, count=int(flag[1] // 0.01))
            self.p.time_step()

if __name__ == '__main__':
    # Инициализация
    v_ = Variables()
    o = Objects(v=v_)
    o.p.show_rate = 10

    # Интегрирование
    o.integrate(t=o.v.TIME, animate=False)

    # Вывод результатов
    tmp = np.array([np.linalg.norm(o.f.line_difference[0][i]) for i in range(len(o.f.line_difference[0]))])
    print(f"Математическое ожидание: {tmp.mean()}, Среднее отклонение: {tmp.std()}")
    talk_decision(cnd=o.v.IF_TALK)
    plot_all(o)
    # plot_distance(o)  # Можно не комментировать
    # plot_sigmas(o)
