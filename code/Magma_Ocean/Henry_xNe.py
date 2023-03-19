import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker, cm, colors
from functions_const import *

SMALL_SIZE = 13
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}


N = 1000
T_surf = np.linspace(1000, 6000, N)
P_surf = np.logspace(-11, 5, N) #MPa
P_bar = P_surf*10

X, Y = np.meshgrid(T_surf, P_surf)
Z = -lnk_Ne_Earth(X, Y)[0]

X_2, Y_2 = np.meshgrid(T_surf, P_bar)

fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
origin = 'lower'


lev_exp = np.arange(1, 60, 1)
levs = np.log10(lev_exp)
cs = ax.contourf(Y_2, X, Z, 100, vmin = 5, vmax = 60, cmap = 'RdYlGn')
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


def x_g(T, P):
    S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
    f_g = P*x_Ne
    return f_g/S_Ne


Z_x = x_g(X,Y)
print(Z_x)

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)

min_xne = -25
max_xne = -20

print('max', np.floor(np.log10(Z_x.min())-1))

print(Z_x.max(),Z_x.min())

lev_exp = np.arange(np.floor(np.log10(Z_x.min())-1),
                    np.ceil(np.log10(Z_x.max())+1), 0.2)
levs = np.power(10, lev_exp)
css = ax.contourf(Y_2, X, Z_x, levs, norm=matplotlib.colors.LogNorm(), cmap = 'magma')
#ee = np.arange(np.floor(np.log10(Z_x.min())-1), np.ceil(np.log10(Z_x.max())+1) , 2)
#levss = np.power(10., ee)

#applyplotstyle(ax)
'''
CS4 = ax.contour(Y_2, X, Z_x, levss,
                 colors=('black',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Magma Ocean concentration',fontsize=title_size)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=contour_label_font_size)
'''
levss = np.arange(np.floor(np.log10(Z_x.min())-1),
                    np.ceil(np.log10(Z_x.max())+1), 2)
levsv = np.power(10, levss)
cbar = fig.colorbar(css, ticks=levsv, norm=colors.LogNorm())
cbar.ax.set_yscale('log')
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Concentration', rotation=270)


ax.set_xlabel('Total Surface Pressure (Bar)')
ax.set_ylabel('Surface Temperature (K)')
ax.set_xscale('log')
ax.set_title(r'Magma Ocean Concentration ($x_{\rm Ne}$) ',fontsize=16)


#ax.xaxis.set_major_locator(mticker.LogLocator(numticks=999))
#ax.xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
plt.tick_params(direction='in', length=6, width=1.5)
plt.tick_params(axis="y", which='major', direction="in", width=1.5, length=8, pad=4)
plt.tick_params(axis="y", which='minor', direction="in", width=1.5, length=4, pad=4)
plt.tick_params(axis="x", which='major', direction="in", width=1.5, length=8, pad=7)
plt.tick_params(axis="x", which='minor', direction="in", width=1.5, length=4, pad=7)
plt.tick_params(right=True, top=True, labelright=False, labeltop=False)
plt.savefig(f"Magma Ocean Concentration", dpi=300, bbox_inches='tight')
plt.show()
