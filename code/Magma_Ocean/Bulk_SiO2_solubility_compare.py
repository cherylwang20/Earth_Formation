import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker, cm, colors

SMALL_SIZE = 15
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}

P_fix = 0.1/10 #MPa

alpha_Ne = 0.401
beta_Ne = -43.36
nu_ca = 13.30
nu_ca_MgO = 7.57
nu_s = -7.18
nu_s_MgO = -2.2
lambda_sio2 = -770.13
lambda_MgO = 2118.9
k_sio2 = -2.06e-6
k_MgO = 4.7e-5
V_Ne = 26.5
Mol_SiO2 = 60.09
Mol_MgO = 40.3
T_ref = 1300  # degree
P_ref = 0.1  # MPa
x_Ne = 2.1e-5  # from Olsen paper #molar ratio


def lnk_Ne_Earth(T, P):
    mol_sio2 = 0.446
    mol_MgO = 0.554
    fun = nu_ca*mol_sio2 + nu_ca_MgO*mol_MgO + nu_s + (lambda_sio2* mol_sio2 + mol_MgO*lambda_MgO) * (1 / T - 1 / T_ref) + (k_sio2*mol_sio2 + k_MgO*mol_MgO) * (P - P_ref) * 10
    result = alpha_Ne * (100 - 100 / (49.14 * 112.33 / 271.85) * fun) + beta_Ne
    return result, fun

def lnk_Ne(T, P):
    fun = nu_ca + nu_s + lambda_sio2 * (1 / T - 1 / T_ref) + k_sio2 * (P - P_ref) * 10
    result = alpha_Ne * (100 - 100 / (Mol_SiO2 * 107.67 / 240.36) * fun) + beta_Ne
    return result, fun

def Ne_res(M_core, P, T, lnk):
    # Ne_con = 3.42e-14
    S_Ne = np.exp(-lnk(T, P)[0])
    #print(S_Ne)
    return M_core * P / S_Ne * x_Ne * 5.972e+24 * 1000  # g

T_list = np.linspace(1500, 4000, 100)

Ne_Pure = Ne_res(1, P_fix, T_list, lnk_Ne)
Ne_Bulk = Ne_res(1, P_fix, T_list, lnk_Ne_Earth)

fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
plt.plot(T_list, Ne_Pure, label = 'Solubility of Pure SiO2')
plt.plot(T_list, Ne_Bulk, label = 'Solubility of Bulk Silicate Earth')
plt.legend()
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Solubility of Neon (g)')
plt.title('Solubility of two Magma Ocean composition under P = 0.1 bar Surface Pressure')
plt.savefig(f"Solubility of two composition at P = {P_fix*10}Bar.png", dpi=300, bbox_inches='tight')
plt.show()