"""НЕ УДАЛЯТЬ! Кое-кто может лишиться парочки конечностей!"""
import scipy
import numpy as np
import matplotlib.pyplot as plt

Radius_orbit = 6800e3
mu = 5.972e24 * 6.67408e-11
w = np.sqrt(mu / Radius_orbit ** 3)
v_orb = np.sqrt(mu / Radius_orbit)
h_orb = Radius_orbit - 6371e3

# Задание параметров
random_rate = 0.01
post_factor = 5
show_factor = 5
r_noise = 0.0
n = 10
dt = 1.
t = [60.*10*(i+0) for i in range(5)]  # ВОЗМОЖНО ТУТ НАДО С НУЛЯ
colors = ["violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
          "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
          "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
          "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
          "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
          "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
          "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
          "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine",
          "violet", "deepskyblue", "maroon", "gold", "forestgreen", "indigo", "olivedrab",
          "slategray", "pink", "salmon", "tan", "steelblue", "peru", "aquamarine"]

# Генерация движения
def get_rand_rvs(r_spread: float = 100, v_spread: float = 0.01):
    r_ = np.random.uniform(-r_spread, r_spread, 3)
    v_ = np.random.uniform(-v_spread, v_spread, 3)
    c_ = np.random.uniform(0, 0.05**2, 1)
    return np.append(np.append(r_, v_), c_).tolist()
def get_atm_params(h: float):
    T = -131.21 + 0.00299 * (h + h_orb)
    p = 2.488 * ((T + 273.1) / 216.6) ** -11.388
    rho = p / (0.2869 * (T + 273.1))
    return rho, T, p
def get_orb_acceleration(r, v):
    return np.array([-2 * w * v[2],
                     -w ** 2 * r[1],
                     2 * w * v[0] + 3 * w ** 2 * r[2]])
def get_full_acceleration(c_resist: float, square: float, r: any, v: any, m: float):
    # force = - r * self.mu / np.linalg.norm(r)
    rho = get_atm_params(r[2])[0]
    force = get_orb_acceleration(r, v)
    v_real = v + np.array([v_orb, 0, 0])
    force -= v_real * np.linalg.norm(v_real) * c_resist * rho * square / 2 / m
    return force
def get_integrate(rvs: any):
    r_ = np.array(rvs[0:3])
    v_ = np.array(rvs[3:6])
    c_ = rvs[6]
    f = get_full_acceleration(c_resist=1.17, square=c_, m=0.03, r=r_, v=v_)
    v_ += f * dt
    r_ += v_ * dt
    return np.append(np.append(r_, v_), c_).tolist()
rvs = []
rvs_rand = []
for i in range(n):
    rvs += get_rand_rvs()
    rvs_rand += get_rand_rvs()

# Расчёт измерений
R = []
rvs_copy = rvs.copy()
for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
    r = []
    time = j * dt
    for i in range(n):
        rvs_copy[i*7:(i+1)*7] = get_integrate(rvs_copy[i*7:(i+1)*7])
    for t1 in t:
        if time == t1:
            for i in range(n):
                r += rvs_copy[i*7:i*7+3]
            for i1 in range(n):
                for i2 in range(i1):
                    dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                    R += [dist_real * (1 + np.random.uniform(-r_noise, r_noise))]


# Отображение
x, y, z = ([], [], [])
rvs_copy = rvs.copy()
ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
for i in range(n):
    x_, y_, z_ = ([], [], [])
    for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
        time = j * dt
        rvs_copy[i*7:(i+1)*7] = get_integrate(rvs_copy[i*7:(i+1)*7])
        x_ += [rvs_copy[i * 7 + 0]]
        y_ += [rvs_copy[i * 7 + 1]]
        z_ += [rvs_copy[i * 7 + 2]]
        for t1 in t:
            if time == t1:
                x += [rvs_copy[i*7+0]]
                y += [rvs_copy[i*7+1]]
                z += [rvs_copy[i*7+2]]
    ax.plot(x_, y_, z_, c=colors[i])
ax.scatter(x, y, z)
plt.show()

def local_func(rvs_: any):
    loss = 0
    counter = 0
    rvs_copy = rvs_.copy().tolist()
    for j in np.linspace(0, t[-1], int(t[-1] // dt), endpoint=False):
        r = []
        time = j * dt
        for i in range(n):
            rvs_copy[i * 7:(i + 1) * 7] = get_integrate(rvs_copy[i * 7:(i + 1) * 7])
        for t1 in t:
            if time == t1:
                for i in range(n):
                    r += rvs_copy[i*7:i*7+3]
                for i1 in range(n):
                    for i2 in range(i1):
                        R_ = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                        loss += (R_ - R[counter]) ** 2
                        counter += 1
    return loss / counter


x0 = np.array(rvs_rand) * random_rate + np.array(rvs) * (1 - random_rate)
res = scipy.optimize.minimize(local_func, x0, tol=1e-4, method='trust-constr',
                              options={'verbose': 3})
print(f"Результатики: {res.x}\nА надо бы {rvs}\n")

# Повторное отображение
ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
x, y, z = ([], [], [])
rvs_copy = rvs.copy()
rvs_res = res.x.copy()
a = [[] for i in range(n)]
for i in range(n):
    x_, y_, z_ = ([], [], [])
    x_1, y_1, z_1 = ([], [], [])
    for j in np.linspace(0, t[-1]*post_factor, show_factor * int(t[-1]*post_factor // dt), endpoint=False):
        time = j * dt
        rvs_res[i*7:(i+1)*7] = get_integrate(rvs_res[i*7:(i+1)*7])
        x_ += [rvs_res[i * 7 + 0]]
        y_ += [rvs_res[i * 7 + 1]]
        z_ += [rvs_res[i * 7 + 2]]
        rvs_copy[i*7:(i+1)*7] = get_integrate(rvs_copy[i*7:(i+1)*7])
        x_1 += [rvs_copy[i * 7 + 0]]
        y_1 += [rvs_copy[i * 7 + 1]]
        z_1 += [rvs_copy[i * 7 + 2]]
        a[i] += [np.linalg.norm(np.array(rvs_res[i*7:i*7+3]) - np.array(rvs_copy[i*7:i*7+3]))]
        for t1 in t:
            if time == t1:
                x += [rvs_copy[i*7+0]]
                y += [rvs_copy[i*7+1]]
                z += [rvs_copy[i*7+2]]
    ax.plot(x_, y_, z_, c=colors[i])
    ax.plot(x_1, y_1, z_1, c=colors[i], ls=":")
ax.scatter(x, y, z)
ax.set_title("Траектории дейтвительные и расчитанные")
ax.set_xlabel("x, м")
ax.set_ylabel("y, м")
ax.set_zlabel("z, м")
plt.show()

for i in range(n):
    plt.plot(np.linspace(0, t[-1]*post_factor, show_factor * int(t[-1]*post_factor // dt), endpoint=False),
             a[i], c=colors[i])
plt.title("Реальная ошибка")
plt.xlabel("Время, с")
plt.ylabel("Ошибка, м")
plt.show()

# Расчёт и отображение измерений
loss = []
counter = 0
rvs_res = res.copy()
for j in np.linspace(0, t[-1]*post_factor, int(t[-1]*post_factor // dt), endpoint=False):
    tmp = 0.
    r = []
    time = j * dt
    for i in range(n):
        rvs_res[i * 7:(i + 1) * 7] = get_integrate(rvs_res[i * 7:(i + 1) * 7])
    for t1 in t[0:-1]:
        if time == t1:
            for i in range(n):
                r += rvs_res[i*7:i*7+3]
            for i1 in range(n):
                for i2 in range(i1):
                    # print(f"{counter+1} / {len(R)}, t = {j}")
                    R_ = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                    tmp += (R_ - R[counter]) ** 2
                    counter += 1
    loss += [np.sqrt(tmp) / counter]
for j in t:
    plt.plot([j, j], [np.min(loss), np.max(loss)], c='gray')
plt.plot([t[0], t[-1]*post_factor], [0, 0], c='gray')
plt.plot(np.linspace(0, t[-1]*post_factor, int(t[-1]*post_factor // dt), endpoint=False), loss)
plt.title("Расчётная ошибка")
plt.xlabel("Время, с")
plt.ylabel("Ошибка, м")
plt.show()
