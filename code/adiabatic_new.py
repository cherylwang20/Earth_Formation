import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

T_1 = 1500 #K
gamma = 9./7. #1.424
R_gas = 8.314 * 10**7
N = 100

def T_adi(a, M, T_top, r, gamma):
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm

    R_out = min(R_H, R_B)
    r_min = R_E * (M/M_earth) ** 0.25


    B = (gamma - 1 )/ gamma * mu / R_gas * G * M
    T = T_disk_chiang(a) - ( B/r - B/r_min )

    return T

# for a 0.1 Earth mass obkect
T_array = [0]*N
r = np.linspace(0.5**0.25*R_E, G * 0.5*M_earth / Cs_disk(T_disk_chiang(1)) ** 2, N)
for i in range(N):
    T_array[i] = T_adi(1, 0.5 * M_earth, T_disk_chiang(1), r[i], gamma)
    if T_array[i] < T_disk_chiang(1):
        T_array[i] = T_disk_chiang(1)

def rho_adi(T_profile):
    S_Ne = np.exp(-lnk_Ne_Earth(T_0, 0.1)[0])
    P_Ne = Ne_reserve * S_Ne / x_Ne / 0.5 / M_earth * 1e6
    rho_0 = P_Ne / Cs_disk(T_profile(1))
    return rho_0*( np.asarray(T_profile)/T_0)**(1/(gamma - 1)), rho_0
#
fig, ax1 = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)

plt.plot((r_max(0.5,T_disk_chiang(1)) - r)/R_E, T_array)
plt.xscale('log')
plt.xlabel('Radius (R_E)')
plt.ylabel('Temperature (K)')
#plt.yscale('log')
plt.show()

rho_neb = 6*10**(-6)

# plt.plot(r/R_E, rho_adi(T_array)[0]/rho_neb)
# plt.xlabel('Atmospheric Radius (R_E)')
# plt.ylabel(r'Density ($\rho/\rho_{neb}$)')
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

def M_adi(x, T_0, M, a, T_top):
    S_Ne = np.exp(-lnk_Ne_Earth(T_0, 0.1)[0])
    #print(f'the mass is {M/M_earth}')
    P_Ne = Ne_reserve * S_Ne / x_Ne / (M/M_earth) / M_earth * 1e6
    #print(P_Ne)
    rho_0 = P_Ne / Cs_disk(T_0)**2
    #print(rho_0)

    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm
    R_out = min(R_H, R_B)
    r_min = R_E * (M / M_earth) ** 0.25

    B = (gamma - 1) / gamma * G * M / Cs_disk(T_top) ** 2
    fun = (1 + B*(1/r_min - 1/R_out))**(1/(gamma-1))
    rho_out = rho_0/fun

    T = T_0 + B / r - B / r_min
    return 4 * np.pi * rho_out * x ** 2 * (1 + B*(1/x - 1/R_out))**(1/(gamma - 1))


M_co = np.logspace(-1, 0, N)
print(f'Mass array:{M_co}')
GCR_adi = [0]*N

def r_min(M):
    return R_E * M ** 0.25

for i in range(N):
    GCR_adi[i] = quad(M_adi, r_min(M_co[i]), r_max(M_co[i], T_1), args=(T_1, M_co[i]*M_earth, 1, T_disk_chiang(1) ))[0]/M_co[i]/M_earth


#print(r_min(0.5), r_max(0.5, T_0))

#print(quad(M_adi, r_min(0.5), r_max(0.5, T_0), args=(T_0, 0.5 * M_earth, 1, T_disk_chiang(1) )))
T_range = np.linspace(1500, 2500, N)



A, B = np.meshgrid(T_range, M_co)
C = np.zeros((N, N))
for j in range(N):
    for i in range(N):
        C[j][i] = quad(M_adi, r_min(M_co[j]), r_max(M_co[j], T_range[i]), args=(T_range[i], M_co[j]*M_earth, 1, T_disk_chiang(1) ))[0]/M_co[j]/M_earth



max_GCR = 4
min_GCR = -10

fig , ax = plt.subplots(figsize=(13, 11), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_GCR, max_GCR, 0.2)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(A, B, C, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Surface Temperature (K)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
#ax.set_xscale('log')
ee = np.arange(min_GCR, max_GCR, 0.2)
levss = np.power(10., ee)

CS4 = ax.contour(A, B, C, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Adiabatic Atmospheric Profile for Chiang Disk Temperature at SMA = 1AU, $\gamma$ = {gamma:.3f}',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"Adiabatic_Chiang_GCR sma = 1AU, 97", dpi=300, bbox_inches='tight')
plt.show()

plt.plot(M_co, GCR_adi)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('GCR')
plt.xlabel('Core Mass (M_Earth)')
plt.show()

