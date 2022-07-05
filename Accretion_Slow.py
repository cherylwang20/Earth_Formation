from functions_const import *

slow_t = 3 * 10**10 #year
yeartos = 3.154e+7

N = 100
sma = np.logspace(-1, 1, N)
Mp = np.logspace(-1, 0, N)

def M_slow(M_core, a):
    return ((M_core*M_earth)**2*M_earth)**(1/3)*(slow_t*yeartos*(a/10)**3)**(-1)

X, Y = np.meshgrid(Mp, sma)
Z = M_slow(X,Y)

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)

max_rate = 18
min_rate = 8

origin = 'lower'
lev_exp = np.arange(min_rate, max_rate, 0.2)
levs = np.power(10, lev_exp)
cs = ax.contourf(Y, X, Z, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Semi-major axis(AU)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_rate, max_rate, 1)
levss = np.power(10., ee)

CS4 = ax.contour(Y, X, Z, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f"Accretion Rate of planetesmial under slow regime",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(r'Accretion Rate', rotation=270)
plt.savefig(f"Mass Accretion Rate_Slow", dpi=300, bbox_inches='tight')
plt.show()


def lum_slow(M, a):
    Rp = R_E * M ** 0.25
    return G*M*M_earth*M_slow(M,a)/Rp

X, Y = np.meshgrid(Mp, sma)
Z = lum_fun(X,Y, M_slow)
fig , ax = plt.subplots(figsize=(10, 8), dpi=80)

max_lum = 28
min_lum = 20
min_lum = 20

origin = 'lower'
lev_exp = np.arange(min_lum, max_lum, 0.2)
levs = np.power(10, lev_exp)
cs = ax.contourf(Y, X, Z, levs, norm=colors.LogNorm(), cmap='Spectral')
ax.set_xlabel('Semi-major axis(AU)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_lum, max_lum, 1)
levss = np.power(10., ee)

CS4 = ax.contour(Y, X, Z, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f"Luminosity of M_core under slow regime",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(r'Luminosity', rotation=270)
plt.savefig(f"Luminosity_Slow", dpi=300, bbox_inches='tight')
plt.show()

def T_surf_slow(M, a):
    T_0 = T_disk_chiang(a)
    L_acc = lum_slow(M,a)
    R_B = G * M *M_earth  / Cs_disk(T_0) ** 2
    return (L_acc/(16*np.pi*(R_B*0.001)**2*sigma))**(1/4)

X, Y = np.meshgrid(Mp, sma)
Z = T_surf_slow(X,Y)

print(Z)

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)

max_T = 5
min_T = 0

origin = 'lower'
lev_exp = np.arange(min_T, max_T, 0.2)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(Y, X, Z, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Semi-major axis(AU)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_T, max_T, 1)
levss = np.power(10., ee)

CS4 = ax.contour(Y, X, Z, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f"T_surf under slow regime",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(r'T_surf', rotation=270)
plt.savefig(f"T_surf_Slow", dpi=300, bbox_inches='tight')
plt.show()