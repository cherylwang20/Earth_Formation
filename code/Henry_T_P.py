import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker, cm, colors
from functions_const import *

SMALL_SIZE = 13
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}

def lnk_Ne(T, P):
    fun = nu_ca + nu_s + lambda_sio2 * (1 / T - 1 / T_ref) + k_sio2 * (P - P_ref) * 10
    result = alpha_Ne * (100 - 100 / (Mol_SiO2 * 107.67 / 240.36) * fun) + beta_Ne
    return result, fun

print(lnk_Ne(5000, 10**(-10)))

T_surf = np.linspace(1000, 6000, 1000)
P_surf = np.logspace(-11, 5, 1000) #MPa
P_bar = P_surf*10

X, Y = np.meshgrid(T_surf, P_surf)
Z = -lnk_Ne(X, Y)[0]
print(Z)
print(np.log10(np.exp(Z)))

X_2, Y_2 = np.meshgrid(T_surf, P_bar)

fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
origin = 'lower'

    # lev_exp = np.arange(np.floor(np.log10(Z.min())-1), np.ceil(np.log10(Z.max())+1))
lev_exp = np.arange(12, 13, 0.1)
levs = np.log10(lev_exp)
cs = ax.contourf(Y_2, X, Z, 100, vmin = 11, vmax = 13, cmap = 'RdYlGn')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ee = np.arange(12, 13, 0.2)
levss = np.power(10., ee)

'''
CS4 = ax.contour(Y_2, X, np.exp(Z), ee,
                     colors=('red',),
                     linewidths=(1,),
                     origin=origin)
'''
ax.set_xlabel('Total Surface Pressure (Bar)')
ax.set_ylabel('Surface Temperature (K)')
ax.set_xscale('log')
ax.set_title('Solubility Constant ',fontsize=16)
#ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)
cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(r'$-\ln(k_{\rm Ne})$', rotation=270)
plt.xticks(size=14)
plt.yticks(size=14)
plt.tick_params(direction='in', length=6, width=1.5)
plt.tick_params(axis="y", which='major', direction="in", width=1.5, length=8, pad=4)
plt.tick_params(axis="y", which='minor', direction="in", width=1.5, length=4, pad=4)
plt.tick_params(axis="x", which='major', direction="in", width=1.5, length=8, pad=7)
plt.tick_params(axis="x", which='minor', direction="in", width=1.5, length=4, pad=7)
plt.tick_params(right=True, top=True, labelright=False, labeltop=False)
plt.savefig(f"Solubility Constant", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

# plt.plot(P_test,-lnk_Ne_Earth(1700, P_test)[0])
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('Pressure (MPa)')
# plt.ylabel('-ln(Ne)')
# #plt.savefig(f"Pressure-solubility constant", dpi=300, bbox_inches='tight')
# #plt.show()
#
# plt.plot(P_test, k_sio2 * (P_test - P_ref) * 10)
# plt.xscale('log')
# #plt.yscale('log')
# plt.xlabel('Pressure (MPa)')
# plt.ylabel('k_SiO2 (P - Pref)')
# #plt.savefig(f"Pressure-k_SiO2", dpi=300, bbox_inches='tight')
# #plt.show()

# actual Ne_reservoir: 64.5e15, we want to see the core_mass and pressure to dissolve other amount of Ne.

# from Iacono-Marziano paper #the henry constant
x_Ne = 2.1e-5  # from Olsen paper #molar ratio

N = 600
P_new = np.logspace(-4, 1, N) #in bar
P_Mpa_new = P_new / 10  #in MPa
M_co = np.logspace(-1, 0.3, N) * 0.7

M_Ne = 20.1797


# total_mole = np.asarray(Ne_reserve)/M_Ne/x_Ne #mol
# total_mole_mass = x_H2*total_mole*M_H2 + x_He*total_mole*M_He + x_Ne*total_mole*M_Ne

def N_to_surf(x):
    return x/x_Ne

def surf_to_N(x):
    return x*x_Ne

def Ne_Plot(T):

    def Ne_res(M_core, P):
        # Ne_con = 3.42e-14
        S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
        #print(S_Ne)
        return M_core * P / S_Ne * x_Ne * 5.972e+24 * 1000  # g

    Mean_S_Ne = np.mean(np.exp(-lnk_Ne_Earth(T, P_Mpa_new)[0]))

    X, Y = np.meshgrid(M_co, P_Mpa_new)
    Z = Ne_res(X, Y)

    #print(Z)

    X_2, Y_2 = np.meshgrid(M_co, P_new)

    fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
    origin = 'lower'

    # lev_exp = np.arange(np.floor(np.log10(Z.min())-1), np.ceil(np.log10(Z.max())+1))
    lev_exp = np.arange(12, 20, 0.2)
    levs = np.power(10, lev_exp)
    cs = ax.contourf(Y_2, X, Z, levs, norm=colors.LogNorm())
    # cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
    ax.set_xlabel('Total Surface Pressure (Bar)')
    ax.set_ylabel('Mantle Mass (M_earth)')
    ax.set_yscale('log')
    ax.set_xscale('log')

    ee = np.arange(11, 19, 1)
    levss = np.power(10., ee)

    CS4 = ax.contour(Y_2, X, Z, levss,
                     colors=('red',),
                     linewidths=(1,),
                     origin=origin)

    Ne_con = 3.42e-14
    Ne_total = 2.735 * 10 ** (15)  # g #Ne_con*M_co*M_Ne*5.972e+27
    # Ne_pre =Ne_total/M_co*k_Ne_1500/x_Ne/5.972e27*1e-5
    Ne_ATM = Ne_con * M_Ne * 5.972e+27  # g
    Earth_Pres = Ne_ATM * Mean_S_Ne / x_Ne / 5.972e27 * 10
    # print(f'{Ne_ATM:.3e}')
    print(f'{Earth_Pres:.3e}')

    plt.plot(Earth_Pres, 0.7, marker="o", markersize=10, markeredgecolor="red", markerfacecolor="aqua")

    # manual_locations = [(-1, -1.4), (-0.62, -0.7), (-2, 0.5), (1.7, 1.2), (2.0, 1.4), (2.4, 1.7)]

    ax.set_title(f'Dissolvable Neon in Deplete Mantle at T = {T+ 273.15}K')
    ax.clabel(CS4, fmt='%2.1e', colors='black', fontsize=12)

    secax = ax.secondary_xaxis('top', functions=( surf_to_N,N_to_surf))
    secax.set_xlabel('Neon Surface Pressure (Bar)')

    cbar = fig.colorbar(cs)
    cbar.ax.get_yaxis().labelpad = 25
    cbar.ax.set_ylabel('Dissolvable Ne (g)', rotation=270)
    print(T)
    T_out = T
    plt.savefig(f"Bulk_Contour_Core_mass for Reservoir at T = {T  + 273.15}K.png", dpi=300, bbox_inches='tight')
    plt.show()


T_list = [1226.85, 1726.85, 2226.85, 2726.85, 3226.85]

#for i in range(len(T_list)):
    #Ne_Plot(T_list[i])

print(1/k_sio2)