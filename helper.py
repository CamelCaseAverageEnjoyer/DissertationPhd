"""НЕ УДАЛЯТЬ! Функция тестирует солвер ЦЭшек по измерениям"""
import scipy
import numpy as np
import matplotlib.pyplot as plt

Radius_orbit = 6800e3
mu = 5.972e24 * 6.67408e-11
w = np.sqrt(mu / Radius_orbit ** 3)

# Задание параметров
random_rate = 0.5
post_factor = 5
show_factor = 1000
n = 10
r_noise = 0.0
t = [60*4*(i+1) for i in range(10)]
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
def get_rand_c(r_spread: float = 100, v_spread: float = 0.01):
    x, y, z = (np.random.uniform(-r_spread, r_spread)), (np.random.uniform(-r_spread, r_spread)), \
        (np.random.uniform(-r_spread, r_spread))
    vx, vy, vz = (np.random.uniform(-v_spread, v_spread)), (np.random.uniform(-v_spread, v_spread)), \
        (np.random.uniform(-v_spread, v_spread))
    # return [2*z + vx/w, vz/w, -3*z - 2*vx/w, x - 2*vz/w, vy/w, y]
    return [2*z + vx/w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y]
def r_hkw(C, w, t):
    return [-3 * C[0] * w * t + 2 * C[1] * np.cos(w * t) - 2 * C[2] * np.sin(w * t) + C[3],
            C[5] * np.cos(w * t) + C[4] * np.sin(w * t),
            2 * C[0] + C[2] * np.cos(w * t) + C[1] * np.sin(w * t)]
C = []
for i in range(n):
    C += get_rand_c()

# Расчёт измерений
R = []
for j in t:
    r = []
    for i in range(n):
        r += [np.array(r_hkw(C[i*6:(i+1)*6], w, j))]
    for i1 in range(n):
        for i2 in range(i1):
            dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
            R += [dist_real * (1 + np.random.uniform(-r_noise, r_noise))]


# Отображение
x, y, z = ([], [], [])
ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
for i in range(n):
    for j in t:
        tmp = r_hkw(C[i*6:(i+1)*6], w, j)
        x += [tmp[0]]
        y += [tmp[1]]
        z += [tmp[2]]
    x_, y_, z_ = ([], [], [])
    for j in np.linspace(t[0], t[-1], 100):
        tmp = r_hkw(C[i*6:(i+1)*6], w, j)
        x_ += [tmp[0]]
        y_ += [tmp[1]]
        z_ += [tmp[2]]
    ax.plot(x_, y_, z_, c=colors[i])
ax.scatter(x, y, z)
plt.show()

def local_func(C_):
    loss = 0
    counter = 0
    for j in t:
        r = []
        for i in range(n):
            r += [np.array(r_hkw(C_[i*6:(i+1)*6], w, j))]
        for i1 in range(n):
            for i2 in range(i1):
                R_ = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                loss += (R_ - R[counter])**2
                counter += 1
    return loss / counter


x0 = np.random.uniform(-100, 100, n * 6) * random_rate + np.array(C) * (1 - random_rate)
res = scipy.optimize.minimize(local_func, x0, tol=1e-4, method='trust-constr',
                              options={'verbose': 3})
print(f"Результатики: {res.x}\nА надо бы {C}\n")

# Повторное отображение
ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')
x, y, z = ([], [], [])
for i in range(n):
    for j in t:
        tmp = r_hkw(C[i*6:(i+1)*6], w, j)
        x += [tmp[0]]
        y += [tmp[1]]
        z += [tmp[2]]
    x_, y_, z_ = ([], [], [])
    for j in np.linspace(t[0], t[-1]*post_factor, show_factor):
        tmp = r_hkw(C[i*6:(i+1)*6], w, j)
        x_ += [tmp[0]]
        y_ += [tmp[1]]
        z_ += [tmp[2]]
    ax.plot(x_, y_, z_, c=colors[i])
    x_, y_, z_ = ([], [], [])
    for j in np.linspace(t[0], t[-1]*post_factor, show_factor):
        tmp = r_hkw(res.x[i*6:(i+1)*6], w, j)
        x_ += [tmp[0]]
        y_ += [tmp[1]]
        z_ += [tmp[2]]
    ax.plot(x_, y_, z_, c=colors[i], ls=":")
ax.scatter(x, y, z)
ax.set_title("Траектории дейтвительные и расчитанные")
ax.set_xlabel("x, м")
ax.set_ylabel("y, м")
ax.set_zlabel("z, м")
plt.show()

for i in range(n):
    a = []
    for j in np.linspace(t[0], t[-1]*post_factor, 100):
        tmp = np.array(r_hkw(C[i*6:(i+1)*6], w, j)) - np.array(r_hkw(res.x[i*6:(i+1)*6], w, j))
        a += [np.linalg.norm(tmp)]
    plt.plot(np.linspace(t[0], t[-1]*post_factor, 100), a, c=colors[i])
plt.title("Реальная ошибка")
plt.xlabel("Время, с")
plt.ylabel("Ошибка, м")
plt.show()

# Расчёт и отображение измерений
loss = []
for j in np.linspace(t[0], t[-1]*post_factor, 100):
    counter = 0
    tmp = 0
    r = []
    r_ = []
    for i in range(n):
        r += [np.array(r_hkw(C[i*6:(i+1)*6], w, j))]
        r_ += [np.array(r_hkw(res.x[i*6:(i+1)*6], w, j))]
    for i1 in range(n):
        for i2 in range(i1):
            R = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
            R_ = np.linalg.norm(r_[i1]) if i1 == i2 else np.linalg.norm(r_[i1] - r_[i2])
            tmp += (R_ - R)**2
            counter += 1
    loss += [np.sqrt(tmp) / counter]
for j in t:
    plt.plot([j, j], [np.min(loss), np.max(loss)], c='gray')
plt.plot([t[0], t[-1]*post_factor], [0, 0], c='gray')
plt.plot(np.linspace(t[0], t[-1]*post_factor, 100), loss)
plt.title("Расчётная ошибка")
plt.xlabel("Время, с")
plt.ylabel("Ошибка, м")
plt.show()
