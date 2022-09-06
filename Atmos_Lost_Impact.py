import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

N = 1000
m_ratio = np.linspace(0, 5.5, N)  # m/M

def X_loss(v_r,m):
    b = v_r* m
    loss = 0.4 * b + 1.8 * b ** 2 - 1.2 * b ** 3
    for i in range(len(loss)):
        if loss[i] > 1:
            loss[i:] = 1
    return loss

#for an adiabatic temperature profile
v_ratio = np.linspace(0.2, 3, 7) #v_imp/v_esc
print(v_ratio)

fig , ax = plt.subplots(figsize=(13, 11), dpi=80)
for i in range(7):
    plt.plot(m_ratio, X_loss(v_ratio[i],m_ratio), label = f'v_imp/v_esc = {v_ratio[i]:.2f}')
plt.xscale('log')
plt.ylim([0,1])
plt.xlim([0,5])
plt.xlabel('m/M')
plt.ylabel(r"$\chi$_loss")
plt.legend()
plt.savefig(f"Giant_Impact_Loss", dpi=300, bbox_inches='tight')
plt.show()

gamma = 1.424
rho_neb = 6*10**(-6)

def M_adi(x, T_0, M, a, T_top):
    B = (gamma - 1) / gamma * G * M / Cs_disk(T_top) ** 2
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm
    R_out = min(R_H, R_B)
    r_min = R_E * (M / M_earth) ** 0.25
    fun = rho_neb*(1 + B*(1/x - 1/R_out))**(1/(gamma - 1))
    #print(rho_0)
    return 4 * np.pi * x**2 * fun

n = 5
M_co = np.logspace(-1, 0, n)
#print(f'Mass array:{M_co}')
Mass_adi = [0]*n

def r_min(M):
    return R_E * M ** 0.25

for i in range(n):
    Mass_adi[i] = quad(M_adi, r_min(M_co[i]), r_max(M_co[i], 1500), args=(1500, M_co[i]*M_earth, 1, T_disk_chiang(1) ))[0]


#for v_imp/v_esc = 1

# m_ratio_2 = (0, 1, 8)
# fig , ax = plt.subplots(figsize=(13, 11), dpi=80)
#
# ax.set_yscale('log')
# ax.set_xlim([0,1])
# ax.set_xlabel('m/M')
# ax.set_ylabel("Total Atmospheric Loss (g)")


def Neon_total(x):
    return x*x_Ne


def total_Neon(ax1):
    y1, y2 = ax1.get_ylim()
    ax_twin.set_ylim(Neon_total(y1), Neon_total(y2))
    ax_twin.figure.canvas.draw()


fig, ax1 = plt.subplots(figsize=(13, 11), dpi=80)

ax_twin = ax1.twinx()

ax1.callbacks.connect("ylim_changed", total_Neon)
for i in range(n):
    ax1.plot(m_ratio, X_loss(1,m_ratio)*Mass_adi[i],label = f'Mass Core = {M_co[i]:.2f}')
    ax1.annotate(f"M = {M_co[i]:.2f}M$_\oplus$", [np.median(m_ratio),np.median(X_loss(1,m_ratio)*Mass_adi[i])+1])
ax1.set_xlim(0, 1)

ax1.set_ylabel('Total Atmospheric Loss (g)')
ax1.set_xlabel('m/M')
ax1.set_yscale('log')
ax_twin.set_ylabel('Neon Loss (g)')
ax_twin.set_yscale('log')
#ax1.legend()
plt.savefig(f"Neon_Loss_Giant Impact", dpi=300, bbox_inches='tight')
plt.show()
#axi.set_yticklabels(Neon_label)

