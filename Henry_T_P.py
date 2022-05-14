import numpy as np
import matplotlib.pyplot as plt

alpha_Ne = 0.401
beta_Ne = -43.36
nu_ca = 13.30
nu_s = -7.18
lambda_sio2 = -770.13
k_sio2 = -2.06e-6
V_Ne = 26.5
M_Ne = 60.09
T_ref = 1300  # degree
P_ref = 0.1  # MPa


def lnk_Ne(T, P):
    fun = nu_ca + nu_s + lambda_sio2 * (1 / T - 1 / T_ref) + k_sio2 * (P - P_ref) * 10
    result = alpha_Ne * (100 - 100 / (M_Ne * 107.67 / 240.36) * fun) + beta_Ne
    return result


print(lnk_Ne(1227, 1.02334021e-05))

P_list = np.logspace(-4, 1, 100)
lnk_Ne_list = lnk_Ne(1727, P_list)
plt.loglog(P_list, -lnk_Ne(1727, P_list))
# plt.xscale('log')

plt.show()
