import matplotlib.pyplot as plt
import numpy as np
#from scipy import interpolate
#import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from functions_const import *

def T_shrink(a, T = 4000, r = 2.5, A = 0.0, M_s = 0.5):
    s = r * R_sun * (1/0.5)**0.4 
    a = a*autocm
    T_eq = T*np.sqrt(s/2/a)*(1 - A)**(1/4)*(1/M_s)**(0.18)
    return T_eq



N = 100
W = 100
M_co = np.linspace(0.1, 1, N)
semi = np.linspace(0.1,  1.0, W)
gamma = 7/5
T_s = 1300

#setting up initial conditions

R_bot = np.zeros((W,N))
R_shrink = np.zeros((W,N))
R_ini = np.zeros((W,N))

rho_out = np.zeros((W, N))
rho_neb = np.zeros((W, N))

P_bot = np.zeros((W, N))
P_ini = np.zeros((W, N))

Mass_adi = np.zeros((W, N))

#print(Mass_adi.shape)

#initial condition for magma ocean ingassing
for j in range(N):
    for i in range(W):
        Mass_adi[i][j] = quad(M_adi, r_min(M_co[j]), r_max(M_co[j], T_disk_chiang(semi[i])), args=(gamma, M_co[j]*M_earth, semi[i], T_disk_chiang(semi[i]), MMSN_V(semi[i]) ))[0]
        P_ini[i][j] = adi_P_T(gamma, M_co[j]*M_earth, semi[i], T_disk_chiang(semi[i]), MMSN_V(semi[i]))[1]

rho_out = np.zeros((W,N))


def r_out_shrink(M, a):
    T_eq = T_shrink(a)
    B  = (gamma - 1) / gamma * G * M / (Cs_disk(T_eq) ** 2)
    r_min = R_E * (M / M_earth) ** 0.25
    R_out = 1/(1/r_min - (T_s/T_eq - 1)/B)

    R_H = a * (M  / 3 / M_sun) ** (1 / 3) * autocm # cm
    R_B = G * M  / Cs_disk(T_disk_chiang(a)) ** 2   # cm
    R_ini = min(R_H, R_B)

    
    
    return R_out, R_ini, r_min
    

for j in range(N):
    for i in range(W):
        R_shrink[i][j], R_ini[i][j], R_bot[i][j] = r_out_shrink(M_co[j]*M_earth, semi[i])



def r_int(x, gamma, M,  T_eq, R_out):
    B  = (gamma - 1) / gamma * G * M / (Cs_disk(T_eq) ** 2)
    return 4*np.pi * x **2 * (1 + B*(1/x - 1/R_out))**(1/(gamma - 1))



for i in range(W):
    for j in range(N):
        rho_out[i][j] =  Mass_adi[i][j]/quad(r_int, R_bot[i][j], R_shrink[i][j], args=(gamma, M_co[j]* M_earth, T_shrink(semi[i]), R_shrink[i][j]))[0]
        rho_neb[i][j] = MMSN_V(semi[i])


def rho_in(gamma, M, a, rho_o):
    T_eq = T_shrink(a)
    B  = (gamma - 1) / gamma * G * M / (Cs_disk(T_eq) ** 2)
    r_min = R_E * (M / M_earth) ** 0.25
    R_out = 1/(1/r_min - (T_s/T_eq - 1)/B)
    rho_i  = rho_o*(1 + B*(1/r_min - 1/R_out))**(1/(gamma - 1))
    return rho_i



for i in range(W):
    for j in range(N):
        P_bot[i][j] = rho_in(gamma, M_co[j]*M_earth, semi[i], rho_out[i][j])*Cs_disk(T_s)**2
        #if P_bot[i][j] < P_ini[i][j]:
            #print(f'New P_bot is smaller than P_ini at {semi[i]}AU, {M_co[j]}M_E')

#print(P_bot)
#print(P_ini)



fig , ax = plt.subplots(figsize=(13, 10), dpi=80)

X, Y = np.meshgrid(semi, M_co)
Z = P_bot*1E-6


lev_exp = np.arange(np.floor(np.log10(Z.min())-1),
                    np.ceil(np.log10(Z.max())+1), 0.2)
levs = np.power(10, lev_exp)
css = ax.contourf(Y, X, Z, levs, norm=matplotlib.colors.LogNorm(), cmap = 'YlOrBr')

applyplotstyle2(ax)

ax.set_title(f'Surface Pressure after Cooling (Bar)',fontsize=title_size)
levss = np.arange(np.floor(np.log10(Z.min())-1),
                    np.ceil(np.log10(Z.max())+1), 2)
levsv = np.power(10, levss)
cbar = fig.colorbar(css, ticks = levsv)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_yscale('log')
cbar.ax.set_ylabel('Pressure (Bar)', rotation=270)
ax.set_ylim([0.1 ,1])
plt.savefig(f"Surface Pressure after Cooling", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()

fig , ax = plt.subplots(figsize=(13, 10), dpi=80)

Z_i = P_ini*1E-6


lev_exp = np.arange(np.floor(np.log10(Z.min())-1),
                    np.ceil(np.log10(Z.max())+1), 0.2)
levs = np.power(10, lev_exp)
css = ax.contourf(Y, X, Z_i, levs, norm=matplotlib.colors.LogNorm(), cmap = 'YlOrBr')

applyplotstyle2(ax)

ax.set_title(f'Surface Pressure during Accretion (Bar)',fontsize=title_size)
levss = np.arange(np.floor(np.log10(Z.min())-1),
                    np.ceil(np.log10(Z.max())+1), 2)
levsv = np.power(10, levss)
cbar = fig.colorbar(css, ticks = levsv)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_yscale('log')
cbar.ax.set_ylabel('Pressure (Bar)', rotation=270)
ax.set_ylim([0.1 ,1])
plt.savefig(f"Surface Pressure during Accretion", dpi=300, bbox_inches='tight')
#plt.show()
plt.close()
