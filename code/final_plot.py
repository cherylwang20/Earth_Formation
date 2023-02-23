import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from functions_const import *


SMALL_SIZE = 20
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)


color = "plasma"
line_c = "black"

#calculate the atmospheric mass as a function of core mass and nebular density
N = 100
semi = 0.1 #AU
M_co = np.linspace(0.1, 1, N)
rho = np.logspace(-18, -5, N)
#print(np.log10(MMSN_V(semi)*1e-8))
semi_out = int(semi * 10.)
gamma = 1.424
T_surf = 1500
avge_mol = 28.9647 #g/mol

Mass_adi = np.zeros((N, N))
Y, X = np.meshgrid(M_co, rho/MMSN_V(semi))

P_bot_adi = np.zeros((N, N))
T_bot_adi = np.zeros((N, N))
M_thres = 0

for j in range(N):
    for i in range(N):
        Mass_adi[i][j] = quad(M_adi, r_min(M_co[j]), r_max(M_co[j], T_surf), args=(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i] ))[0]
        P_bot_adi[i][j] = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[1]
        check_t = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[3]
     # check if the temperature if able the melting point
        T_bot_adi[i][j] = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[3]
        if check_t < 1700:
            M_thres = M_co[j]

print(M_thres)
M_thres_list = [M_thres]*len(M_co)

max_atm = 21
min_atm = 6

#now we calculate the total amount of dissolvable gas bassed on this temperature and pressure

def Ne_res(M_core, P, T):
    # Ne_con = 3.42e-14
    S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
    # print(S_Ne)
    return M_core * P / S_Ne * x_Ne * 5.972e+24 * 1000  # g

Neon_dis = np.zeros((N, N))
for j in range(N):
    for i in range(N):
        Neon_dis[i][j] = Ne_res(M_co[j]*0.7,P_bot_adi[i][j]*10**(-7), T_bot_adi[i][j]) #g

#print(Neon_dis)

fig , ax = plt.subplots(figsize=(12, 12), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_atm,max_atm, 0.2)
levs = np.power(10, lev_exp)

cs = ax.contourf(X, Y, Neon_dis, levs, norm=colors.LogNorm(), cmap=color)
ee = np.arange(min_atm, max_atm, 1)
levss = np.power(10., ee)
applyplotstyle(ax)
CS4 = ax.contour(X, Y, Neon_dis, levss,
                 colors=(line_c,),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Neon Capacity (g)',fontsize=title_size)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=contour_label_font_size)

#cbar = fig.colorbar(cs)
#cbar.ax.get_yaxis().labelpad = 25
#cbar.ax.set_ylabel('Neon Capacity(g)', rotation=270)
ax.set_ylim([0.1,1])
plt.fill_between(rho/MMSN_V(semi),M_thres_list, color = 'grey', hatch = r'\\', alpha = 0.75)






def degas_percent(p):
    bulk_con = 12.5*3.42E-16 #mol/g #97% degassing rate
    Mass_silicate = M_co*M_earth*0.7
    degas_c = bulk_con*Mass_silicate*mole_Ne/(1 - p/100)
    return degas_c


total_neon = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        total_neon[i][j] =Neon_dis[i][j]


rho_des_999 = des_rho(total_neon, degas_percent(99.9))
rho_des_945 = des_rho(total_neon,degas_percent(94.5))
#rho_des_75= des_rho(total_neon,degas_percent(75))
# = des_rho(total_neon, degas_percent(75))

fig , ax = plt.subplots(figsize=(12, 12), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_atm,max_atm, 0.2)
levs = np.power(10, lev_exp)

M_embryos = 0.25
#eccentricity = 0.01
ecc = [0, 1E-3, 1E-2]
boun = [np.min(rho)]*len(rho)
cs = ax.contourf(X, Y, total_neon, levs, norm=colors.LogNorm(), cmap=color,rasterized=True)


plt.plot(rho_des_999/MMSN_V(semi),M_co,'w--',linewidth = 5, label = '99.9%')
plt.plot(rho_des_945/MMSN_V(semi),M_co,'c--',linewidth = 5, label = '94.5%')
#plt.plot(rho_des_40/MMSN_V(semi), M_co, 'w--', linewidth = 5, label = '40%')
ax.set_ylim([0.1,1])

rho_lower_limit = [10**(-8)]*len(M_co)
rho_upper_limit = [1]*len(M_co)
plt.fill_betweenx(M_co, rho_des_945/MMSN_V(semi), rho_des_999/MMSN_V(semi),color = "w", alpha = 0.25, hatch=r"//")

plt.fill_betweenx(M_co, rho_upper_limit, rho_des_999/MMSN_V(semi),color = "yellow", alpha = 0.3, hatch=r"*")
plt.fill_betweenx(M_co, rho_des_945/MMSN_V(semi), rho_lower_limit,color = "red", alpha = 0.3, hatch=r"o")


M_p = np.asarray([0.25, 0.5])
ratio = np.asarray([8, 15.75])
plt.fill_between(rho/MMSN_V(semi),M_thres_list, color = 'grey', hatch = r'\\', alpha = 0.9)


for j in range(2):
    for i in range(len(ecc)):
        if i == 0:
            plt.plot(max_rho(semi, M_p[j], ecc[i], tx(M_embryos,M_sun,semi,ecc[i], ratio[j]))/MMSN_V(semi), M_p[j],marker = "o", markersize = 12, markerfacecolor = "blue", markeredgecolor = 'blue')
        if i == 1:
            plt.plot(max_rho(semi, M_p[j], ecc[i], tx(M_embryos,M_sun,semi,ecc[i], ratio[j]))/MMSN_V(semi), M_p[j],marker = "o", markersize = 12, markerfacecolor = "orange", markeredgecolor = 'orange')
        if i == 2:
            plt.plot(max_rho(semi, M_p[j], ecc[i], tx(M_embryos,M_sun,semi,ecc[i], ratio[j]))/MMSN_V(semi), M_p[j],marker = "o", markersize = 12,  markerfacecolor = "green", markeredgecolor = 'green')

legend_elements = [Line2D([0], [0], color = 'w', ls = '-', lw=4, label='99.9%'),
                    Line2D([0], [0], color = 'cyan', ls = '-', lw=4, label='94.5%'),
                    Line2D([0], [0], marker='o', color='blue', label='e = 0',
                          markerfacecolor='blue', markersize=12),
                    Line2D([0], [0], marker='o', color='orange', label='e = 0.001',
                          markerfacecolor='orange', markersize=12),
                    Line2D([0], [0], marker='o', color='green', label='e = 0.01',
                          markerfacecolor='green', markersize=12),
                   ]

def ftotau(x):
    V = -tau_disk*np.log(x)
    return [abs(int(z)) for z in V]

def tautof(x):
    return np.exp(x/tau_disk)
#ax2 = ax.secondary_xaxis("top", functions=(ftotau,tautof))
applyplotstyle(ax)
ax2 = ax.twiny()

ax2.set_xbound(ax.get_xbound())
array = np.logspace(-8,0,9)
print(array)
print(ftotau(array))
ax2.set_xticks(np.linspace(0, 1, 9))
ax2.set_xticklabels(ftotau(array))
ax2.set_xlabel("Dissipation Time (Myr)", fontsize = contour_label_font_size + 3, labelpad=10)
ee = np.arange(min_atm, max_atm , 1)
levss = np.power(10., ee)

CS4 = ax.contour(X, Y, total_neon, levss,
                 colors=(line_c,),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Available Dissolvable Neon (g) a = {semi}AU',fontsize=title_size, y=1.1)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize = contour_label_font_size)
tickss = np.arange(min_atm, max_atm, 1)

ax.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.1, 1.0),
          fancybox=True, fontsize = legend_font_size)
plt.savefig(f"Available Neon under an adiabatic profile, a = {semi_out}AU, D'Alessio", dpi=300, bbox_inches='tight')
plt.show()
