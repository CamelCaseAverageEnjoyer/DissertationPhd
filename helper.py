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

pass
