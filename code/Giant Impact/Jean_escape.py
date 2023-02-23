import matplotlib.pyplot as plt
import numpy as np
from sympy.solvers import solve
from sympy import Symbol
x = Symbol('x')

from functions_const import *

N = 100
M_co = np.linspace(0.1, 1, N)
T = 1500#T_disk(0.5)
n_ex = 3.*10**6


#print(M_co/r_exo**2)

def F_esc(m):
    r_exo = ((G * M_earth * M_co * m / k_B / T / np.sqrt(2) / n_ex / sigma_neon) ** (1 / 2))
    g_r = G * M_co * M_earth / r_exo ** 2
    f_esc = G*M_co*M_earth*m/k_B/T/r_exo
    F_r = 2 / np.sqrt(np.pi) * np.sqrt(2 * k_B * T / m) * (1 + f_esc) * np.exp(-f_esc) * n_ex
    return F_r


A_planet = 4*np.pi*r_min(M_co)**2
#print(F_r*m_Ne*A_planet)


T_range = np.linspace(100, 2500, 100)
v_t_Ne = np.sqrt(k_B*T_range/(m_Ne))
v_t_He = np.sqrt(k_B*T_range/(m_He))
v_t_H = np.sqrt(k_B*T_range/(m_H))

def v_to_M(x):
    return (18*(x)**2*R_E/G/M_earth**(1/4))**(4/3)/M_earth

def M_to_v(x):
    return 1/6*np.sqrt(2*G*(x*M_earth)**(3/4)/R_E/M_earth**(1/4))


##plotting the 1/6 escape speed
fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

ax.plot(T_range, v_t_H, label = 'Hydrogen')
ax.plot(T_range, v_t_He, label = 'Helium')
ax.plot(T_range, v_t_Ne, label = 'Neon')
ax.set_xlabel('Temperature (K)')
ax.set_ylabel('1/6 Escape Velocity (cm/s)')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_ylim([10**4, 7*10**5])
ax.legend()


secax = ax.secondary_yaxis('right', functions=(v_to_M, M_to_v))
secax.set_ylabel('Core Mass (M_Earth)')
plt.savefig(f"Helium & Neon Escape Velocity", dpi=300, bbox_inches='tight')
plt.show()

fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

plt.plot(M_co, F_esc(m_H), label = 'Hydrogen')
plt.plot(M_co, F_esc(m_He), label = 'Helium')
plt.plot(M_co, F_esc(m_Ne), label = 'Neon')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Core Mass (M_Earth)')
plt.ylabel(r'Escape Flux ($cm^{-2}s^{-1}$)')
plt.legend()
plt.savefig(f"Escape Flux", dpi=300, bbox_inches='tight')
plt.show()

def Mass_esc(flux,mole):
    A_planet = 4*np.pi*r_min(M_co)**2
    t_E = 4.6 * 10 ** 9 * 3.154* 10**7
    mass = A_planet*flux*t_E/av*mole
    return mass


fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

plt.plot(M_co, Mass_esc(F_esc(m_H),mole_H), label = 'Hydrogen')
plt.plot(M_co, Mass_esc(F_esc(m_He),mole_He), label = 'Helium')
plt.plot(M_co, Mass_esc(F_esc(m_Ne),mole_Ne), label = 'Neon')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Core Mass (M_Earth)')
plt.ylabel(r'Mass Loss (g)')
plt.legend()
plt.savefig(f"Mass Loss", dpi=300, bbox_inches='tight')
plt.show()
