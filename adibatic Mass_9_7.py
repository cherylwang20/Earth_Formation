import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import matplotlib.ticker as mticker

from functions_const import *

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

for j in range(N):
    for i in range(N):
        Mass_adi[i][j] = quad(M_adi, r_min(M_co[j]), r_max(M_co[j], T_surf), args=(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i] ))[0]
        P_bot_adi[i][j] = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[1]
        T_bot_adi[i][j] = adi_P_T(gamma, M_co[j]*M_earth, semi, T_disk(semi), rho[i])[3]

#desire_rho = np.loadtxt(open("desire_rho.csv"), delimiter=",")
#rho_fun = interpolate.interp1d(M_co[2:],desire_rho[2:],fill_value='extrapolate')
#desire_rh = np.loadtxt(open("desire_rho_9.csv"), delimiter=",")
#rho_fun_9 = interpolate.interp1d(M_co[2:],desire_rh[2:],fill_value='extrapolate')


#rho_des = rho_fun(M_co) #99.9% degassing rate
#rho_des_9 = rho_fun_9(M_co) # 99% degassing rate

fig , ax = plt.subplots(figsize=(12, 10), dpi=80)

max_atm = 27
min_atm = 7


origin = 'lower'
lev_exp = np.arange(min_atm,max_atm, 0.2)
levs = np.power(10, lev_exp)

cs = ax.contourf(X, Y, Mass_adi, levs, norm=colors.LogNorm(), cmap=color)
applyplotstyle(ax)

ee = np.arange(min_atm,max_atm, 1)
levss = np.power(10., ee)

CS4 = ax.contour(X, Y, Mass_adi, levss,
                 colors=(line_c,),
                 linewidths=(1,),
                 origin=origin)


ax.set_title(f'Atmospheric Mass of Adiabatic Profile',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)
cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Atmospheric Mass (g)', rotation=270)
plt.savefig(f"Atmospheric Mass of Adiabatic Profile, {semi_out}AU", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()

#convert the atmospheric mass to the Neon mass

Neon_atm = Mass_adi/avge_mol*x_Ne*mole_Ne

fig , ax = plt.subplots(figsize=(12, 10), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_atm,max_atm, 0.2)
levs = np.power(10, lev_exp)

cs = ax.contourf(X, Y, Neon_atm, levs, norm=colors.LogNorm(), cmap=color)
ee = np.arange(min_atm,max_atm, 1)
levss = np.power(10., ee)

CS4 = ax.contour(X, Y, Neon_atm, levss,
                 colors=(line_c,),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Neon Mass of Adiabatic Profile',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)

applyplotstyle(ax)

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Neon Mass (g)', rotation=270)
plt.savefig(f"Neon Mass of Adiabatic Atmosphere, a = {semi_out}AU", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()

# we calculate the surface pressure and surface temperature using the adiabtic profile
# and nebular density with a certain disk profile
#print(f'The surface temperature is: {T_bot_adi}')


#Now we calculate the surface pressure and temperature of an adiabatic atmosphere
#by using Henry's law, and the inferred adiabtic index, we have

fig , ax = plt.subplots(figsize=(12, 10), dpi=80)

max_press = -4
min_press = -14

lev_exp = np.arange(min_press,max_press, 0.2)
levs = np.power(10, lev_exp)

cs = ax.contourf(X, Y, P_bot_adi*1E-6, levs, norm=colors.LogNorm(), cmap='RdYlGn')
ee = np.arange(min_press, max_press, 1)
levss = np.power(10., ee)

applyplotstyle(ax)

CS4 = ax.contour(X, Y, P_bot_adi*1E-6, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Surface Pressure of Adiabatic Atmosphere',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Surface Pressure (bar)', rotation=270)
plt.savefig(f"Surface Pressure of Adiabatic Atmosphere, a = {semi_out}AU", dpi=300, bbox_inches='tight')
plt.show()
plt.close()


fig , ax = plt.subplots(figsize=(12, 10), dpi=80)

max_temp = 5000
min_temp = 1000

lev_exp = np.arange(min_press,max_press, 0.2)
levs = np.power(10, lev_exp)

cs = ax.contourf(X, Y, T_bot_adi, 100, vmin = 1000, vmax = 5000, cmap='RdYlGn')
ee = np.arange(min_press, max_press, 1)
levss = np.power(10., ee)
applyplotstyle(ax)

ax.set_title(f'Surface Temperature of Adiabatic Atmosphere',fontsize=16)

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Surface Temperature (K)', rotation=270)
plt.savefig(f"Surface Temperature of Adiabatic Atmosphere, a = {semi_out}AU", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

#print(P_bot_adi) #ba g/cm/s^2


#now we calculate the total amount of dissolvable gas bassed on this temperature and pressure

def Ne_res(M_core, P, T):
    # Ne_con = 3.42e-14
    S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
    # print(S_Ne)
    return M_core * P / S_Ne * x_Ne * 5.972e+24 * 1000  # g

Neon_dis = np.zeros((N, N))
for j in range(N):
    for i in range(N):
        Neon_dis[i][j] = Ne_res(M_co[j]*0.7,P_bot_adi[i][j], T_bot_adi[i][j]) #g

#print(Neon_dis)

fig , ax = plt.subplots(figsize=(12, 10), dpi=80)

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
ax.set_title(f'Neon Capacity (g)',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Neon Capacity(g)', rotation=270)
plt.savefig(f"Neon Capacity under an adiabatic profile, a = {semi_out}AU", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()

def degas_percent(p):
    bulk_con = 12.5*3.42E-16 #mol/g #97% degassing rate
    Mass_silicate = M_co*M_earth*0.7
    degas_c = bulk_con*Mass_silicate*mole_Ne/p*100
    return degas_c

# fig , ax = plt.subplots(figsize=(12, 10), dpi=80)
# plt.plot(bulk_con*Mass_silicate*mole_Ne*97/3,M_co, label = '97%')
# plt.plot(bulk_con*Mass_silicate*mole_Ne*3,M_co, label = '75%')
# plt.plot(degas_99,M_co, label = '99.9%')
# plt.xlabel('Total Neon in Mantle (g)')
# plt.ylabel('Core Mass (M_Earth)')
# plt.xlim([10**14,10**18])
#
# plt.legend()
# plt.savefig(f"Theortical 20Ne concentration", dpi=300, bbox_inches='tight')
# #plt.show()
# plt.close()

total_neon = np.zeros((N, N))

for i in range(N):
    for j in range(N):
        total_neon[i][j] = min(Neon_dis[i][j], Neon_atm[i][j])


rho_des_999 = des_rho(total_neon, degas_percent(99.9))
rho_des_99 = des_rho(total_neon,degas_percent(99))
rho_des_75 = des_rho(total_neon, degas_percent(75))

fig , ax = plt.subplots(figsize=(12, 10), dpi=80)
origin = 'lower'
lev_exp = np.arange(min_atm,max_atm, 0.2)
levs = np.power(10, lev_exp)

M_embryos = 0.5
#eccentricity = 0.01
ecc = [0, 1E-3, 1E-2]
boun = [np.min(rho)]*len(rho)
cs = ax.contourf(X, Y, total_neon, levs, norm=colors.LogNorm(), cmap=color,rasterized=True)
#plt.plot([MMSN_V(semi)]*len(M_co),M_co, 'y-',label = 'MMSN', linewidth = 4)
plt.plot(rho_des_999/MMSN_V(semi),M_co,'w--',linewidth = 5, label = '99.9%')
plt.plot(rho_des_99/MMSN_V(semi),M_co,'c--',linewidth = 5, label = '99%')
plt.plot(rho_des_75/MMSN_V(semi), M_co, 'g--', linewidth = 5, label = '75%')

plt.fill_betweenx(M_co, boun, rho_des_999/MMSN_V(semi),color = "w", alpha = 0.25, hatch=r"//")
for i in range(len(ecc)):
    plt.plot(max_rho(semi, M_embryos, ecc[i], tx(M_embryos,M_sun,semi,ecc[i]))/MMSN_V(semi), M_embryos,marker = "o", markersize = 10, label = f'e = {ecc[i]}')#, markerfacecolor = "aqua", markeredgecolor = 'red')
print(max_rho(semi, M_embryos, 0.01, tx(M_embryos,M_sun,semi,0.01)))
applyplotstyle(ax)
ee = np.arange(min_atm, max_atm, 1)
levss = np.power(10., ee)

CS4 = ax.contour(X, Y, total_neon, levss,
                 colors=(line_c,),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Available Dissolvable Neon (g) a = {semi}AU',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)
tickss = np.arange(min_atm, max_atm, 1)
cbar = fig.colorbar(cs)

cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Available Neon(g)', rotation=270)
cbar.ax.tick_params(axis = 'both', which  = 'both', length = 0, labelsize=13)
ax.legend(loc='upper center', bbox_to_anchor=(0.1, 1.0),
          fancybox=True, fontsize = 13)
plt.savefig(f"Available Neon under an adiabatic profile, a = {semi_out}AU, D'Alessio", dpi=300, bbox_inches='tight')
plt.show()
#plt.close()

#np.savetxt('output_9.csv',total_neon,delimiter=",")
#np.savetxt('degas_9.csv',degas_99,delimiter=",")
## now we perform an interpolation between the density and the Neon.
#degas_99 is the amount of target line we want to draw
