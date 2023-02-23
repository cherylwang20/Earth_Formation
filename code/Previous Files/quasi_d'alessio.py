from functions_const import *

N = 100
T_surf = np.linspace(1500, 3000, N)

sma = np.logspace(-1, 0.1, N)
M_co = np.logspace(-1, 0, N)
Z = np.zeros((N, N))

A, B = np.meshgrid(T_surf, M_co)
C = np.zeros((N, N))

for j in range(len(M_co)):
    for i in range(len(T_surf)):
        C[j][i] = GCR_quasi(1, M_co[j],T_disk,T_surf[i])

max_GCR = 4
min_GCR = -10

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)

origin = 'lower'
lev_exp = np.arange(min_GCR, max_GCR, 0.2)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(A, B, C, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Surface Temperature (K)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_GCR, max_GCR, 1)
levss = np.power(10., ee)

CS4 = ax.contour(A, B, C, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f"Quasi-Isothermal Profile for D'alessio Disk Temperature at SMA = 1AU",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


def tick_function(X):
    V = (X/150)**(-7/3)
    return V #["%.3f" % z for z in V]

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"Quasi-Isothermal_D'Alessio_GCR for sma and M at sma = 1AU", dpi=300, bbox_inches='tight')
plt.show()