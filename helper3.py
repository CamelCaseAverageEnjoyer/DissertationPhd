import matplotlib.pyplot as plt
import numpy as np

POINTS = 33
LINKS = 46
CIRCLE = 0.4
h = 2.
x = [0, h, h, 2*h,  h,  0, -h,   -h,    0,    h,  2*h,  3*h,  4*h,  4*h, 4*h,  5*h,  5*h, -2*h,   -h,    0,    0,    0,    h,    h,  2*h,  3*h,  4*h,  6*h, 6*h, 6*h,  5*h, 6*h, 6*h]
y = [0, h, 0,   0, -h, -h, -h, -2*h, -2*h, -2*h, -2*h, -2*h, -2*h, -1*h,   0, -1*h, -2*h, -3*h, -3*h, -3*h, -4*h, -5*h, -4*h, -3*h, -3*h, -3*h, -3*h, -4*h, 3*h, 4*h, -4*h,  -h, 2*h]

link = [[0, 0, 0.1, True] for _ in range(LINKS)]  # Начало, конец, длина, индикатор
f = open("vla.txt", "r", encoding="utf-8")
counter = 0
flag = 0
for line in f:
    # print(f"line:{line} -> {line.split()}")
    tmp = line.split()[-1]  # "\ufeff1")[0]
    tmp = tmp.replace(',', '.')
    # print(f"tmp={tmp}, counter={counter}, flag={flag}")
    tmp = float(tmp) if flag >= 2 else int(tmp)
    link[counter][flag] = tmp
    counter += 1
    if counter == LINKS:
        counter = 0
        flag += 1
f.close()

f = open("vladvladvlad.txt", "r", encoding="utf-8")
counter = 0
for line in f:
    # print(f"line:{line} -> {line.split()}")
    tmp = line.split()[-1]  # "\ufeff1")[0]
    # print(f"tmp={tmp}, counter={counter}, flag={flag}")
    link[counter][3] = bool(int(tmp))
    counter += 1
f.close()


fig, ax = plt.subplots(figsize=(10, 10))

circles = [plt.Circle((x[i], y[i]), CIRCLE, color='c') for i in range(POINTS)]
for i in range(LINKS):
    c = 'red' if link[i][3] else 'gray'
    ax.plot([x[link[i][0]-1], x[link[i][1]-1]], [y[link[i][0]-1], y[link[i][1]-1]], c=c)
    ax.text((x[link[i][0]-1] + x[link[i][1]-1]) / 2 - 0.4, (y[link[i][0]-1] + y[link[i][1]-1])/2, f'{link[i][2]}')
for i in range(POINTS):
    ax.add_patch(circles[i])
    tmp = 0 if (i+1) < 10 else CIRCLE/3
    ax.text(x[i]-CIRCLE/2 - tmp, y[i]-CIRCLE/2, f'{i+1}')

plt.xlim((-7.5, 17.5))
plt.ylim((-12.5, 12.5))

plt.show()
