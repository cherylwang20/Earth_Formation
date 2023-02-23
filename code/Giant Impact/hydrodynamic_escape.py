import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

N = 1000
#T_homo = 1500
M = np.linspace(0.1, 1, N)
T_i = 10000

b = 2.7*10**17 * T_i**0.75#1.8 * 10**19 #/cm/s
print(b)
t_1 = 1.6
t_2 = 0.1
X_H = 0.912
X_Ne = 0.000035

m_1 = 2 #amu
m_2 = 20.1797

def mc(m, T, t):
    phi = 23.3 * t**(-1.23) * 0.15*2/4
    r = r_min(m)*2
    F = phi*r/G/m/M_earth/m_H/2 #1.5*10**11#
    #print(F)
    fun = k_B*T/b/X_H*F*r**2/G/m/M_earth/m_H
    #print(fun)
    return m_1 + fun

mc_t1 = mc(M, T_i, t_1)
mc_t2 = mc(M, T_i, t_2)

fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

plt.plot(M, mc_t1, label =f'{t_1} Gyr')
plt.plot(M, mc_t2, label =f'{t_2} Gyr')
plt.xscale('log')
plt.xlabel('Core Mass (M_E)')
plt.ylabel('Crossover Mass (amu)')
plt.yscale('log')
plt.legend()
plt.savefig(f"Crossover Mass", dpi=300, bbox_inches='tight')
plt.show()


def F_Ne(m,t):
    phi = 23.3 * t ** (-1.23) * 0.15 * 2 / 4
    r = r_min(m) * 2
    mcross = mc(m, T_i, t)
    #print(f'mass is: {mcross}')
    F = phi * r / G / m / M_earth / m_H / 2  # 1.5*10**11#
    return X_Ne/X_H*F*(mcross - m_2)/(mcross - m_1)

t_int = np.linspace(0.5, 4.6, N)
F_total = [0]*len(M)

for i in range(len(M)):
    F_total[i] = np.sum(F_Ne(M[i],t_int)) #per atom/cm^2/s
    if F_total[i] < 0:
        F_total[i] = 0
    else:
        F_total[i] = F_total[i]/av #mole/cm^2/s
    #print(F_total[i])

fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)
plt.plot(M, F_total)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Core Mass (M_E)')
plt.ylabel('Total Escape Flux (mole/cm^2/s)')

plt.show()


fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

F_Ne_t1 = F_Ne(M,t_1)
F_Ne_t2 = F_Ne(M,t_2)

plt.plot(M, F_Ne_t1, label =f'{t_1} Gyr')
plt.plot(M, F_Ne_t2, label =f'{t_2} Gyr')
plt.xscale('log')
plt.xlabel('Core Mass (M_E)')
plt.ylabel(r'Escape Flux (cm$^{-2}$s$^{-1}$)')
#plt.yscale('log')
plt.legend()
plt.savefig(f"Escape Flux_hydrodynamics", dpi=300, bbox_inches='tight')
plt.show()

def Mass_esc(flux,m, M_core):
    A_planet = 4*np.pi*r_min(M_core)**2
    t_E = (10**8 - 10*7) * 3.154* 10**7
    mass = A_planet*flux*t_E*m
    GCR = mass/M_core/M_earth/0.00058
    return mass, GCR


fig, ax1 = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

print(Mass_esc(F_total,20.2,M)[0])

ax2 = ax1.twinx()
ax1.plot(M, Mass_esc(F_total,20.2,M)[0], 'b--')
ax2.plot(M, Mass_esc(F_total,20.2,M)[1], 'g-')
ax1.set_xscale('log')
ax1.set_xlabel('Core Mass (M_E)')
ax1.set_ylabel(r'Mass Loss of Neon(g)', color = 'b')
ax2.set_ylabel('Gas to Core Ratio (GCR)', color = 'g')
ax2.set_yscale('log')
ax1.set_yscale('log')
plt.savefig(f"Mass Loss & GCR via hydrodynmical Loss", dpi=300, bbox_inches='tight')
plt.show()