# from main_objects import *
# from plot_func import *

# Определение конфигурации согласно статье
'''a_m_f = np.zeros((self.f.n, self.f.n))
a_p_f = np.zeros((self.f.n, self.f.n))
a_m_c = np.zeros((self.f.n, self.c.n))
a_p_c = np.zeros((self.f.n, self.c.n))

def local_func(x1):
    anw = 0
    for i__c in range(self.c.n):
        for i__f in range(self.f.n):
            anw += a_m_c[i__f][i__c] + a_p_c[i__f][i__c]
    for i_1f in range(self.c.n):
        for i_2f in range(self.f.n):
            anw += a_m_f[i_1f][i_2f] + a_p_f[i_1f][i_2f]
    return anw

def local_constraints(x1):
    anw = []
    x = x1.reshape((self.f.n, 3))
    z = np.vstack([np.hstack((np.eye(3), x.T)), np.hstack((x, x @ x.T))])

    for i__c in range(self.c.n):
        for i__f in range(self.f.n):
            e = np.hstack([[int(j == i__f) for j in range(self.f.n)], self.c.r_orf[i__c]])
            anw += [-self.c.calc_dist[i__c][i__f][len(self.c.calc_dist[i__c][i__f]) - 1] ** 2 +
                    e.T @ z @ e - a_p_c[i__f][i__c] + a_m_c[i__f][i__c]]
    for i_1f in range(self.c.n):
        for i_2f in range(self.f.n):
            if i_1f != i_2f:
                e = np.hstack([[int(j == i_1f) - int(j == i_2f) for j in range(self.f.n)], np.zeros(3)])
                anw += [-self.f.calc_dist[i_1f][i_2f][len(self.f.calc_dist[i_1f][i_2f]) - 1] ** 2 +
                        e.T @ z @ e - a_p_f[i_1f][i_2f] + a_m_f[i_1f][i_2f]]
    return anw

tol = 1
opt = {'verbose': 3, 'gtol': 1e-20, 'xtol': 1e-20}
x0 = flatten(self.r_orf_estimation)
res = scipy.optimize.minimize(local_func, x0, tol=tol, method='trust-constr', options=opt,
                              constraints={'type': 'eq', 'fun': local_constraints})
self.r_orf_estimation -= 0'''

import numpy as np
import random
from sympy import *
from pysr import pysr, best, best_callable
Radius_orbit = 6800e3
mu = 5.972e24 * 6.67408e-11
w_hkw = np.sqrt(mu / (Radius_orbit ** 3))
w_0_hkw = w_hkw

Rp_x, Rp_z, Vp_x, Vp_z, wp_y, m, M, J_y, Jp_y = var('Rp_x Rp_z Vp_x Vp_z wp_y m M J_y Jp_y', real=True, constant=True)
x_0, z_0, x_c, z_c, xp_c, zp_c, w_0, phi_0 = var('x_0 z_0 x_c z_c xp_c zp_c w_0 phi_0', real=True, constant=True)
Rc_x, Rc_z, R_x_0, R_z_0, r_x_0, r_z_0, r1_x, r1_z, T = var('Rc_x Rc_z R_x_0 R_z_0 r_x_0 r_z_0 r1_x r1_z T', real=True, constant=True)

v0_x, v0_z, t, w_y = var('v0_x v0_z t w_y', real=True)
R_x, R_z, V_x, V_z, r_x, r_z, v_x, v_z, phi = symbols("R_x R_z V_x V_z r_x r_z v_x v_z phi", cls=Function)

'''v_x_0 = v0_x + Vp_x + wp_y * sin(phi_0)*(x_0 - xp_c)
v_z_0 = v0_z + Vp_z - wp_y * cos(phi_0)*(x_0 - xp_c)
V_x_0 = -v0_x*m/M + Vp_x + wp_y * sin(phi_0)*(x_c - xp_c)
V_z_0 = -v0_z*m/M + Vp_z - wp_y * cos(phi_0)*(x_c - xp_c)
w_y = ((m+M)*(Rp_z*Vp_x - Rp_x*Vp_z) + Jp_y * wp_y - M*(R_z_0*V_x_0 - R_x_0*V_z_0) - m*(r_z_0*v_x_0 - r_x_0*v_z_0))/J_y'''

v_x_0 = v0_x*cos(phi_0) - v0_z*sin(phi_0) + Vp_x + wp_y * (sin(phi_0)*(x_0 - xp_c) + cos(phi_0)*(z_0 - zp_c))
v_z_0 = v0_x*sin(phi_0) + v0_z*cos(phi_0) + Vp_z - wp_y * (cos(phi_0)*(x_0 - xp_c) - sin(phi_0)*(z_0 - zp_c))
V_x_0 = -v0_x*m/M*cos(phi_0) + v0_z*m/M*sin(phi_0) + Vp_x + wp_y * (sin(phi_0)*(x_c - xp_c) + cos(phi_0)*(z_c - zp_c))
V_z_0 = -v0_x*m/M*sin(phi_0) - v0_z*m/M*cos(phi_0) - Vp_z + wp_y * (cos(phi_0)*(x_c - xp_c) - sin(phi_0)*(z_c - zp_c))
w_y = ((m+M)*(Rp_z*Vp_x - Rp_x*Vp_z) + Jp_y * wp_y - M*(R_z_0*V_x_0 - R_x_0*V_z_0) - m*(r_z_0*v_x_0 - r_x_0*v_z_0))/J_y

# eq1 = Eq(phi(t), phi_0 + t * w_y)
eq1 = Eq(phi(t).diff(t), w_y)
eq2 = Eq(R_x(t).diff(t), V_x(t))
eq3 = Eq(R_z(t).diff(t), V_z(t))
eq4 = Eq(V_x(t).diff(t), -2*w_0*V_z(t))
eq5 = Eq(V_z(t).diff(t), 2*w_0*V_x(t) + 3*w_0**2*R_z(t))
eq6 = Eq(r_x(t).diff(t), v_x(t))
eq7 = Eq(r_z(t).diff(t), v_z(t))
eq8 = Eq(v_x(t).diff(t), -2*w_0*v_z(t))
eq9 = Eq(v_z(t).diff(t), 2*w_0*v_x(t) + 3*w_0**2*r_z(t))
# eq10 = Eq(r1_x,  cos(phi(t))*(r_x(T) - Rc_x) + sin(phi(t))*(r_z(T) - Rc_z))

a1, a2, a3, a4, a5, a6, c1, c2, b1, b2, b3 = var('a1 a2 a3 a4 a5 a6 c1 c2 b1 b2 b3', real=True)
A1 = a1 * v0_x + a2 * v0_z + a3
A2 = a4 * v0_x + a5 * v0_z + a6
B = b1 * v0_x + b2 * v0_z + b3
dr_x = A1*cos(B) - A2*sin(B) + c1
dr_z = A1*sin(B) + A2*cos(B) + c2

a1_p = 3
a2_p = 5
a3_p = 7
a4_p = 11
a5_p = 13
a6_p = 17

b1_p = 19
b2_p = 23
b3_p = 29

c1_p = 31
c2_p = 37
dr_x_p = dr_x.subs([(a1, a1_p), (a2, a2_p), (a3, a3_p), (a4, a4_p), (a5, a5_p), (a6, a6_p),
                    (b1, b1_p), (b2, b2_p), (b3, b3_p), (c1, c1_p), (c2, c2_p)])
dr_z_p = dr_z.subs([(a1, a1_p), (a2, a2_p), (a3, a3_p), (a4, a4_p), (a5, a5_p), (a6, a6_p),
                    (b1, b1_p), (b2, b2_p), (b3, b3_p), (c1, c1_p), (c2, c2_p)])
dr_x_p_func = lambdify([v0_x, v0_z], dr_x_p)
dr_z_p_func = lambdify([v0_x, v0_z], dr_z_p)

N_array = 10
x_array = [[random.uniform(-0.05, 0.05), random.uniform(-0.05, 0.05)] for _ in range(N_array)]
y_array = [np.sqrt(dr_x_p_func(x_array[i][0], x_array[i][1])**2 + dr_z_p_func(x_array[i][0], x_array[i][1])**2) for i in range(N_array)]

eq_pysr = pysr(
    x_array,
    y_array,
    niterations=2,
    binary_operators=["+", '-'],
    unary_operators=[],
)
print(eq_pysr)
