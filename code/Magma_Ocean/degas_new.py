from functions_const import *
import matplotlib.pyplot as plt
import numpy as np

T_s = 1200
N = 10
gamma = 1.4
M_co = np.linspace(0.1, 1, N)
a = np.linspace(0.1, 1, N)
rho = np.logspace(-18, -5, N)
semi = 1.0
semi_out = int(semi * 10.)

f_degas = np.zeros((N, N))
f_degas_2 = np.zeros((N,N))
P_s = 10**(-8)

def Ne_res(M_core, P, T):
    # Ne_con = 3.42e-14
    S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
    # print(S_Ne)
    return M_core * P / S_Ne * x_Ne *5.972e+24 * 1000  # g

for j in range(N):
    for i in range(N):
        P_i = adi_P_T(gamma, M_co[j]*M_earth, a[i], T_disk(a[i]), rho[i])[1]*1E-6 #MPa
        T_i = adi_P_T(gamma, M_co[j]*M_earth, a[i], T_disk(a[i]), rho[i])[3]
        P_s = P_bot(gamma,M_co[j]*M_earth, a[i],  T_disk(a[i]), T_s, rho[i])
        N1 = Ne_res(M_co[j], P_i, T_i)
        N2 = Ne_res(M_co[j], P_s, T_s)
        f_degas[i][j] = (N1 - N2)/N1
        #at a fixed distance
        P_j = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[1]*1E-6 #MPa
        T_j = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[3]
        P_ss = P_bot(gamma,M_co[j]*M_earth, semi,  T_disk(semi), T_s, rho[i])
        N11 = Ne_res(M_co[j], P_j, T_j)
        N22 = Ne_res(M_co[j], P_ss, T_s)
        f_degas_2[i][j] = (N11 - N22)/N11

print(f_degas.max(),f_degas.min())
print(f_degas_2.max(),f_degas_2.min())

P_ii = adi_P_T(gamma, 0.5*M_earth, semi, T_disk(semi), 10**(-18))[1]*1E-6 #MPa
T_ii = adi_P_T(gamma, 0.5*M_earth, semi, T_disk(semi), 10**(-18))[3]
Ni = Ne_res(0.5, P_ii, T_ii)
Nj = Ne_res(0.5, P_ii, T_s)
print((Ni - Nj)/Ni)

fig , ax = plt.subplots(figsize=(13, 10), dpi=80)

degas_min = 0
degas_max = 1.05

lev_exp = np.arange(degas_min,degas_max, 0.05)
lev_label = np.arange(degas_min, degas_max, 0.1)
Y, X = np.meshgrid(M_co, a)

cs = ax.contourf(X, Y, f_degas, lev_exp, cmap='spring')

applyplotstyle2(ax)

ax.set_title(f'Degassing Rate of Magma Ocean(%)',fontsize=title_size)

cbar = fig.colorbar(cs, ticks = lev_label )
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Degassing rate (%)', rotation=270)
ax.set_ylim([0.1, 1])
plt.savefig(f"Degassing Rate of Magma Ocean", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

fig , ax = plt.subplots(figsize=(13, 10), dpi=80)

degas_min = 0
degas_max = 1.05

lev_exp = np.arange(degas_min,degas_max, 0.05)
lev_label = np.arange(degas_min, degas_max, 0.1)
Y, X = np.meshgrid(M_co, rho/MMSN_V(semi))

cs = ax.contourf(X, Y, f_degas_2, lev_exp, cmap='summer_r')

applyplotstyle(ax)

ax.set_title(f'Degassing Rate of Magma Ocean at a = {semi}AU (%)',fontsize=title_size)

cbar = fig.colorbar(cs, ticks = lev_label )
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Degassing rate (%)', rotation=270)
ax.set_ylim([0.1, 1])
plt.savefig(f"Degassing Rate of Magma Ocean at a = {semi_out}AU", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

