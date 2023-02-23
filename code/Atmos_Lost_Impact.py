import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

m_ratio = np.linspace(0.1, 2, 1000)  # m/M

def X_loss(v_r,m):
    b = v_r* m
    loss = 0.4 * b + 1.8 * b ** 2 - 1.2 * b ** 3
    for i in range(len(loss)):
        if loss[i] > 1:
            loss[i:] = 1
    return loss

n = 5
#for an adiabatic temperature profile
v_orbit = np.linspace(10,30,n)*2
v_escape = np.sqrt(2*G*0.5*M_earth/r_min(0.5))/100000

v_ratio = v_orbit/v_escape
print(v_ratio)

fig , ax = plt.subplots(figsize=(13, 11), dpi=80)
for i in range(n):
    plt.plot(m_ratio, X_loss(v_ratio[i],m_ratio), label = f'v_imp/v_esc = {v_ratio[i]:.2f}')
#plt.xscale('log')
plt.ylim([0,1])
plt.xlim([0.1,1])
plt.xlabel('m/M')
plt.ylabel(r"$\chi$_loss")
plt.legend()
plt.savefig(f"Giant_Impact_Loss", dpi=300, bbox_inches='tight')
plt.show()

gamma = 1.424
rho_neb = 6*10**(-6)

n = 5
M_co = np.linspace(0.2, 1, n)
#print(f'Mass array:{M_co}')
N = 100
semi = 1.0 #AU
rho = MMSN_V(semi)*1E-7
gamma = 1.424
T_surf = 1500
avge_mol = 28.9647 #g/mol

Mass_adi = [0]*n



for i in range(n):
    Mass_adi[i] = quad(M_adi, r_min(M_co[i]), r_max(M_co[i], T_surf), args=(gamma, M_co[i]*M_earth, semi, T_disk(semi), rho ))[0]

Mass_main = quad(M_adi, r_min(0.5), r_max(0.5, T_surf), args=(gamma, 0.5*M_earth, semi, T_disk(semi), rho ))[0]
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
    ax1.plot(m_ratio, X_loss(v_ratio[i], m_ratio) * Mass_main, label=f'v_impact = {v_orbit[i]:.2f}')
    #ax1.annotate(f"v_impact = {v_orbit[i]:.2f}km/s",[np.median(m_ratio)-0.1, np.median(X_loss(v_ratio[i], m_ratio) * Mass_main) + 1])
ax1.set_xlim(0.1, 2)

ax1.set_ylabel('Total Atmospheric Loss (g)')
ax1.set_xlabel('m/M')
#ax1.set_yscale('log')
ax_twin.set_ylabel('Neon Loss (g)')
#ax_twin.set_yscale('log')
ax1.set_xscale('log')
ax1.legend()
plt.savefig(f"Neon_Loss_Giant Impact", dpi=300, bbox_inches='tight')
plt.show()
#axi.set_yticklabels(Neon_label)

