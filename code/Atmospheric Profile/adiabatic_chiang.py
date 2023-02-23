import matplotlib.pyplot as plt
import numpy as np

from functions_const import *
from scipy.interpolate import interp1d


T_0 = 2000
gamma = 1.424
R_gas = 8.314 * 10**7
N = 1000

gamma_array = np.loadtxt("gamma.txt")
r_array = np.linspace(0.5**0.25*R_E,24.3*0.5**0.25*R_E, len(gamma_array))
gamma_inter = interp1d(r_array, gamma_array, kind='cubic', fill_value='interpolate')

r_array =  np.linspace(0.5**0.25*R_E,24.3*0.5**0.25*R_E, N)

print(24.3*0.5**0.25*R_E/37672511968)
def M_adi(x, a, M, P_bot, T, T_top):
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm

    R_out = min(R_H, R_B)
    r_min = R_E * (M) ** 0.25
    R_out_E = R_out / R_E

    rho_bot = P_bot / Cs_disk(T) ** 2
    B = (gamma - 1/gamma)*mu/R_gas*G*M

    func_new = (T_0 + B/x - B/r_min)**(1/(gamma - 1))* x**2

    C = T_0 - (gamma - 1) / gamma * mu / R_gas * G * M / r_min
    func = x**2 * ((gamma - 1)/gamma * mu/R_gas*G*M/x + C)**(1/(gamma - 1))
    return 4*np.pi*rho_bot*func_new

def T_adi(a, M, T_top, r, gamma):
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm

    R_out = min(R_H, R_B)
    r_min = R_E * (M/M_earth) ** 0.25


    B = (gamma - 1 )/ gamma * mu / R_gas * G * M
    T = T_0 + B/r - B/r_min

    return T

# for a 0.1 Earth mass obkect
r = np.linspace(0.5**0.25*R_E, 37672511968, N)
T_gamma_inter = [0]*len(r_array)


for i in range(len(r_array)):
    T_gamma_inter[i] = T_adi(1, 0.5 * M_earth, T_disk_chiang(1), r_array[i], gamma_inter(r_array)[i])
    if T_gamma_inter[i] < T_disk_chiang(1):
        T_gamma_inter[i] = T_disk_chiang(1)

plt.plot(r_array/R_E, T_gamma_inter)
plt.xscale('log')
#plt.yscale('log')
plt.show()

S_Ne = np.exp(-lnk_Ne_Earth(T_0, 0.1)[0])
P_Ne = Ne_reserve * S_Ne / x_Ne / 0.5 / M_earth * 1e6
rho_0 = P_Ne/Cs_disk(1500)
plt.plot(r_array/R_E, rho_0*( np.asarray(T_gamma_inter)/T_0)**(1/(1.424 - 1)))
plt.xscale('log')
plt.yscale('log')
plt.show()

def GCR_adi(a, M, T_fun, T_bot):
    T_top = T_fun(a)
    r_min = R_E * (M) ** 0.25
    S_Ne = np.exp(-lnk_Ne_Earth(T_bot, 0.1)[0])
    P_Ne = Ne_reserve * S_Ne / x_Ne / M / M_earth * 1e6
    res = quad(M_adi, r_min, r_max(M, T_bot), args=(a, M * M_earth, P_Ne, T_bot, T_top))[0] #Mass given in grams
    return res/M_earth #/ M / M_earth


M_co =  np.logspace(-1, 0, N)

GCR_adiabatic = [0]*len(M_co)
for i in range(len(M_co)):
    GCR_adiabatic[i] = GCR_adi(1, M_co[i],T_disk_chiang,T_0)
plt.plot(M_co, GCR_adiabatic)
plt.xscale('log')
plt.yscale('log')
plt.show()