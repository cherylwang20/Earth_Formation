from functions_const import *

N = 100
T_surf = np.linspace(1500, 3000, N)
T_0 = 1500

sma = np.logspace(-1, 0.1, N)
M_co = np.logspace(-1, 0, N)
Z = np.zeros((N, N))


A, B = np.meshgrid(T_surf, M_co)
C = np.zeros((N, N))

for j in range(len(M_co)):
    for i in range(len(T_surf)):
        C[j][i] = GCR_quasi(1, M_co[j],T_disk_chiang,T_surf[i])

X, Y = np.meshgrid(sma, M_co)

for j in range(len(M_co)):
    for i in range(len(sma)):
        Z[j][i] = GCR_quasi(sma[i], M_co[j] ,T_disk_chiang, 1500)

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)


max_GCR = 4
min_GCR = -10

origin = 'lower'
lev_exp = np.arange(min_GCR, max_GCR, 0.2)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(X, Y, Z, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Semi-major axis(AU)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_GCR, max_GCR, 1)
levss = np.power(10., ee)

CS4 = ax.contour(X, Y, Z, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Quasi-Isothermal Profile for Chiang Disk Temperature at T = {T_0}K',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


def tick_function(X):
    V = (X/150)**(-7/3)
    return V #["%.3f" % z for z in V]



# secax = ax.secondary_xaxis('top', functions = (T_disk_chiang, tick_function))
# secax.set_xlabel('Disk Temperature (K)')

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"Quasi-Isothermal_Chiang_GCR for sma and M at T = {T_0}", dpi=300, bbox_inches='tight')
plt.show()

#T_surf against core mass
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
ax.set_title(f'Quasi-Isothermal Profile for Chiang Disk Temperature at SMA = 1AU',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


def tick_function(X):
    V = (X/150)**(-7/3)
    return V #["%.3f" % z for z in V]

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"Quasi-Isothermal_Chiang_GCR for T_surf and M at sma = 1AU", dpi=300, bbox_inches='tight')
plt.show()