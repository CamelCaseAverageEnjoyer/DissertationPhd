"""НЕ УДАЛЯТЬ! Перебор параметров С по измерениям!"""
import matplotlib.pyplot as plt
from tiny_functions import *
from frases import *

Radius_orbit = 6800e3
mu = 5.972e24 * 6.67408e-11
w = np.sqrt(mu / Radius_orbit ** 3)

# Задание параметров
N = 1000
if_quaternion = False
random_rate = 0.1
post_factor = 1
show_alpha_factor = 1000
show_factor = 100
n = 1
r_noise = 0.0
t = [160*15*(i+1) for i in range(6)]
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
def get_gain(r: any):
    r1 = r / np.linalg.norm(r)
    return 2 - np.linalg.norm(np.dot(r1, np.array([1, 0, 0]))) ** 2 - \
        np.linalg.norm(np.dot(r1, np.array([0, 1, 0]))) ** 2
def get_rand_c(r_spread: float = 100, v_spread: float = 0.01):
    x, y, z = np.random.uniform(-r_spread, r_spread, 3)
    vx, vy, vz = np.random.uniform(-v_spread, v_spread, 3)
    a, b, g = np.random.uniform(-100, 100, 3)
    if if_quaternion:
        return [2*z + vx/w, vz / w, -3 * z - 2 * vx / w, x - 2 * vz / w, vy / w, y, a, b, g]
    return [2*z + vx/w, vz/w, -3*z - 2*vx/w, x - 2*vz/w, vy/w, y]
def r_hkw(C, w, t):
    return [-3 * C[0] * w * t + 2 * C[1] * np.cos(w * t) - 2 * C[2] * np.sin(w * t) + C[3],
            C[5] * np.cos(w * t) + C[4] * np.sin(w * t),
            2 * C[0] + C[2] * np.cos(w * t) + C[1] * np.sin(w * t)]
len_x = 9 if if_quaternion else 6
ev = np.array([1., 0., 0.])
C = []
for i in range(n):
    C += get_rand_c()

# Расчёт измерений
R = []
for j in t:
    r = []
    for i in range(n):
        r += [np.array(r_hkw(C[i*len_x:i*len_x+6], w, j))]
    for i1 in range(n):
        for i2 in range(i1 + 1):
            dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
            if if_quaternion:
                g1 = get_gain(euler2rot_matrix(C[i1*9+6] / 30, C[i1*9+7] / 30, C[i1*9+8] / 30) @ ev)
                g2 = get_gain(euler2rot_matrix(C[i2*9+6] / 30, C[i2*9+7] / 30, C[i2*9+8] / 30) @ ev)
                signal_rate = g1 if i1 == i2 else g1 * g2
                dist_calc = dist_real * np.sqrt(signal_rate)
            else:
                dist_calc = dist_real
            R += [dist_calc * (1 + np.random.uniform(-r_noise, r_noise))]

C_all = [[] for _ in range(N)]
lines = [[] for _ in range(N)]
losses = []
for i_n in range(N):
    print(Fore.LIGHTCYAN_EX + f"{i_n+1}/{N}: {int((i_n+1)/N * 100)}%")
    tmp = []
    for i in range(n):
        tmp += get_rand_c()
    C_all[i_n] = (np.array(tmp) * random_rate + np.array(C) * (1 - random_rate)).tolist()

    loss = 0.
    counter = 0
    for j in t:
        r = []
        for i in range(n):
            r += [np.array(r_hkw(C_all[i_n][i * len_x:i * len_x + 6], w, j))]
        for i1 in range(n):
            for i2 in range(i1 + 1):
                dist_real = np.linalg.norm(r[i1]) if i1 == i2 else np.linalg.norm(r[i1] - r[i2])
                if if_quaternion:
                    g1 = get_gain(
                        euler2rot_matrix(C_all[i_n][i1 * 9 + 6] / 30,
                                         C_all[i_n][i1 * 9 + 7] / 30,
                                         C_all[i_n][i1 * 9 + 8] / 30) @ ev)
                    g2 = get_gain(
                        euler2rot_matrix(C_all[i_n][i2 * 9 + 6] / 30,
                                         C_all[i_n][i2 * 9 + 7] / 30,
                                         C_all[i_n][i2 * 9 + 8] / 30) @ ev)
                    signal_rate = g1 if i1 == i2 else g1 * g2
                    dist_calc = dist_real * np.sqrt(signal_rate)
                else:
                    dist_calc = dist_real
                R_ = dist_calc * (1 + np.random.uniform(-r_noise, r_noise))

                loss += (R_ - R[counter]) ** 2
                counter += 1
    losses += [loss]

# Отображение
ax = plt.figure(figsize=(7, 7)).add_subplot(projection='3d')

for i_n in range(N):
    for i in range(n):
        x_, y_, z_ = ([], [], [])
        for j in np.linspace(t[0], t[-1], show_factor):
            tmp = r_hkw(C_all[i_n][i * len_x:i * len_x + 6], w, j)
            x_ += [tmp[0]]
            y_ += [tmp[1]]
            z_ += [tmp[2]]
        tmp = ((np.max(losses) - losses[i_n]) / np.max(losses))
        if tmp > 0.995:
            a = clip(tmp*show_alpha_factor - show_alpha_factor + 1, 0, 1)
            # ax.plot(x_, y_, z_, c=((1 - a), a, 0, a))
            ax.plot(x_, y_, z_, c=colors[i*2], alpha=a)

x, y, z = ([], [], [])
for i in range(n):
    for j in t:
        tmp = r_hkw(C[i*len_x:i*len_x+6], w, j)
        x += [tmp[0]]
        y += [tmp[1]]
        z += [tmp[2]]
    x_, y_, z_ = ([], [], [])
    for j in np.linspace(t[0], t[-1], show_factor):
        tmp = r_hkw(C[i*len_x:i*len_x+6], w, j)
        x_ += [tmp[0]]
        y_ += [tmp[1]]
        z_ += [tmp[2]]
    ax.plot(x_, y_, z_, c=colors[i*2+1])
ax.scatter(x, y, z)

plt.show()
