import numpy as np

from functions_const import *

T_1 = 1500
gamma = 9./7. #1.424
R_gas = 8.314 * 10**7
N = 100

def T_adi(a, M, T_top, r, gamma):
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm

    R_out = min(R_H, R_B)
    r_min = R_E * (M/M_earth) ** 0.25


    B = (gamma - 1 )/ gamma * mu / R_gas * G * M
    T = T_1 + B/r - B/r_min

    return T

# for a 0.1 Earth mass obkect
T_array = [0]*N
r = np.linspace(0.5**0.25*R_E, G * 0.5*M_earth / Cs_disk(T_disk_chiang(1)) ** 2, N)
for i in range(N):
    T_array[i] = T_adi(1, 0.5 * M_earth, T_disk_chiang(1), r[i], gamma)
    if T_array[i] < T_disk_chiang(1):
        T_array[i] = T_disk_chiang(1)


rho_neb = 6*10**(-6)


def adi_profile(r, T_0, M, a, T_top):
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
    rho_fun = rho_out*(1 + B*(1/r - 1/R_out))**(1/(gamma-1))

    C = (gamma - 1 )/ gamma * mu / R_gas * G * M

    T = T_0 + C / r - C / r_min

    for i in range(len(r)):
        if rho_fun[i] > rho_neb:
            print(i)
            return r[i]/R_E
    for i in range(len(T)):
        if T[i] < T_disk_chiang(1):
            print(i)
            return r[i]/R_E
    return R_out/R_E

T_range = np.linspace(1500, 2500, N)
M_co = np.logspace(-1, 0, N)

r = np.linspace(0.5**0.25*R_E, G * 0.5*M_earth / Cs_disk(T_disk_chiang(1)) ** 2, N)
Q = np.zeros((N, N))

for j in range(N):
    for i in range(N):
        Q[j][i] = adi_profile(r, T_range[i], M_co[j] * M_earth, 1, T_disk_chiang(1))
print(Q)
W, E = np.meshgrid(T_range, M_co)

min_r = 0
max_r = 2
intv = 0.1

fig , ax = plt.subplots(figsize=(13, 11), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_r, max_r, intv)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(W, E, Q, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Surface Temperature (K)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
#ax.set_xscale('log')
ee = np.arange(min_r, max_r, intv*2)
levss = np.power(10., ee)

CS4 = ax.contour(W, E, Q, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Radius to Reach Nebular Temperature at SMA = 1AU, $\gamma$ = {gamma:.3f}',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Radius (R_E)', rotation=270)
plt.savefig(f"Radius to reach nebular T sma = 1AU, 97", dpi=300, bbox_inches='tight')
plt.show()
#print(K)

def rho_out(T_0, M, a, T_top):
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
    #rho_fun = rho_out*(1 + B*(1/r - 1/R_out))**(1/(gamma-1))
    return rho_out

Q = np.zeros((N, N))

for j in range(N):
    for i in range(N):
        Q[j][i] = rho_out(T_range[i], M_co[j] * M_earth, 1, T_disk_chiang(1))
print(Q)
W, E = np.meshgrid(T_range, M_co)

min_r = -15
max_r = -5
intv = 0.1

fig , ax = plt.subplots(figsize=(13, 11), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_r, max_r, intv)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(W, E, Q, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Surface Temperature (K)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
#ax.set_xscale('log')
ee = np.arange(min_r, max_r, intv*2)
levss = np.power(10., ee)

CS4 = ax.contour(W, E, Q, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Boundary Density at SMA = 1AU, $\gamma$ = {gamma:.3f}',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(f'Density (g/cm$^3$)', rotation=270)
plt.savefig(f"Boundary Density sma = 1AU, 97", dpi=300, bbox_inches='tight')
plt.show()