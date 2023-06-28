import plotly.express as px
import plotly.graph_objs as go
import matplotlib.pyplot as plt

from tiny_functions import *
FEMTO_RATE = 1e2
CUBE_RATE = 1e2

def show_chipsat(o, j, clr, opacity):
    global FEMTO_RATE
    x, y, z = ([], [], [0 for _ in range(4)])
    for x_shift in [-o.f.size * FEMTO_RATE, o.f.size * FEMTO_RATE]:
        for y_shift in [-o.f.size * FEMTO_RATE, o.f.size * FEMTO_RATE]:
            x += [x_shift]
            y += [y_shift]
    A = quart2dcm(o.f.q[j])
    for i in range(4):
        r = A @ np.array([x[i], y[i], z[i]])
        x[i] = r[0] + o.f.r[j][0]
        y[i] = r[1] + o.f.r[j][1]
        z[i] = r[2] + o.f.r[j][2]
    xl, yl, zl = ([], [], [])
    '''for i in range(len(o.line_k[j]) // 3):
        xl += [o.line_k[j][3*i + 0]]
        yl += [o.line_k[j][3*i + 1]]
        zl += [o.line_k[j][3*i + 2]]'''
    return [go.Mesh3d(x=x, y=y, z=z, color=clr, opacity=opacity), go.Scatter3d(x=xl, y=yl, z=zl, mode='lines')]

def show_cubesat(o, j):
    global CUBE_RATE
    total_cubes = 30
    r = [[[] for _ in range(total_cubes)] for _ in range(3)]
    # r[x/y/z][0..5 - yellow, 6..29 - gray][1..4 - sides of square]
    sequence = [[0, 0], [0, 1], [1, 0], [1, 1]]
    for i in range(3):
        for k in range(2):
            r[i][k + i * 2] = [(-1)**(k+1) * o.c.size[i] * CUBE_RATE for _ in range(4)]

    shift = [[-o.c.size[i] * CUBE_RATE, o.c.size[i] * CUBE_RATE] for i in range(3)]
    legs = [o.c.legs_x, o.c.legs_x, o.c.legs_z]
    n_legs = int((total_cubes - 6) // 6)
    bound_legs = [[[((-1)**sequence[s][0] * o.c.size[0] - legs[0]) * CUBE_RATE,
                    ((-1)**sequence[s][0] * o.c.size[0] + legs[0]) * CUBE_RATE],
                   [((-1)**sequence[s][1] * o.c.size[1] - legs[1]) * CUBE_RATE,
                    ((-1)**sequence[s][1] * o.c.size[1] + legs[1]) * CUBE_RATE],
                   [(-o.c.size[2] - legs[2]) * CUBE_RATE, (o.c.size[2] + legs[2]) * CUBE_RATE]] for s in range(n_legs)]
    # bound_cube = [[-o.c.size[0], o.c.size[1]], [-o.c.size[1], o.c.size[1]], [-o.c.size[2], o.c.size[2]]]

    for s in range(n_legs):
        # ind1 = 0 + int(s // 2 < 1)
        # ind2 = 1 + int(s // 2 < 2)
        r[0][(s + 1) * 6 + 0] = [bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0], bound_legs[s][0][0]]
        r[1][(s + 1) * 6 + 0] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 0] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 1] = [bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 1] = [bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 1] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 2] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 2] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][0]]
        r[2][(s + 1) * 6 + 2] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 3] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 3] = [bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 3] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][1], bound_legs[s][2][1]]

        r[0][(s + 1) * 6 + 4] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 4] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 4] = [bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0], bound_legs[s][2][0]]

        r[0][(s + 1) * 6 + 5] = [bound_legs[s][0][0], bound_legs[s][0][1], bound_legs[s][0][0], bound_legs[s][0][1]]
        r[1][(s + 1) * 6 + 5] = [bound_legs[s][1][0], bound_legs[s][1][0], bound_legs[s][1][1], bound_legs[s][1][1]]
        r[2][(s + 1) * 6 + 5] = [bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1], bound_legs[s][2][1]]

    for i in range(3):
        for k in range(2):
            # 0 -> 1, 2
            # 1 -> 0, 2
            # 2 -> 0, 1
            tmp = k + i * 2
            ind1 = 0 + int(i < 1)
            ind2 = 1 + int(i < 2)
            print(f"i={i}. k={k}  | t={tmp}, 1={ind1}, 2={ind2}")
            for m in range(4):
                r[ind1][tmp] += [shift[ind1][sequence[m][0]]]
                r[ind2][tmp] += [shift[ind2][sequence[m][1]]]
    '''for s in range(int((total_cubes - 6) // 6)):
        for i in range(3):
            for k in range(2):
                # 0 -> 1, 2
                # 1 -> 0, 2
                # 2 -> 0, 1
                tmp = k + i * 2
                ind1 = 0 + int(i < 1)
                ind2 = 1 + int(i < 2)
                print(f"i={i}. k={k}  | t={tmp}, 1={ind1}, 2={ind2}")
                for m in range(4):
                    # shift[i][sequence[s][0]] + (-1)**(k+1) * legs[i]
                    r[ind1][(s + 1) * 6 + tmp] += [shift[ind1][sequence[s][0]] + sequence[m][0] * legs[i]]
                    r[ind2][(s + 1) * 6 + tmp] += [shift[ind2][sequence[s][1]] + sequence[m][1] * legs[i]]'''
    print(f"Итого r={r}")

    A = quart2dcm(o.c.q[j])
    for k in range(total_cubes):
        for i in range(4):
            # print(f"i={i}. k={k}, r={len(r[0])}")
            r1 = A @ np.array([r[0][k][i], r[1][k][i], r[2][k][i]])
            r[0][k][i] = r1[0] + o.c.r[j][0]
            r[1][k][i] = r1[1] + o.c.r[j][1]
            r[2][k][i] = r1[2] + o.c.r[j][2]
    anw = []
    for i in range(total_cubes):  # range(total_cubes)  [6, 7, 8, 9, 10, 11]
        print(f"x={r[0][i]}")
        print(f"y={r[1][i]}")
        print(f"z={r[2][i]}\n")
        color = 'yellow' if i < 6 else 'gray'
        anw += [go.Mesh3d(x=r[0][i], y=r[1][i], z=r[2][i], color=color, opacity=1)]
        # anw += [go.Mesh3d(x=x[i], y=y[i], z=z[i], color='yellow', opacity=1)]
    print(f"____{legs[0]/o.c.size[0]}")
    return anw

def show_chipsats(o, clr: str = 'lightpink', opacity: float = 1):
    data = []
    for i in range(o.f.n):
        data += show_chipsat(o, i, clr, opacity)
    for i in range(o.c.n):
        data += show_cubesat(o, i)
    return go.Figure(data=data, layout=go.Layout(autosize=False, width=900, height=700))

def plot_all(o):
    f = show_chipsats(o)
    f.show()
