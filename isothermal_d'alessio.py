from functions_const import *

N = 100
sma = np.logspace(-1, 0.1, N)
M_co = np.logspace(-1, 0, N)
Z = np.zeros((N, N))

Y, X = np.meshgrid(M_co, sma)
for j in range(len(M_co)):
    for i in range(len(sma)):
        Z[j][i] = GCR_isothermal(sma[i], M_co[j],T_disk)

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)


max_GCR = 4
min_GCR = -8

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
ax.set_title(f"Isothermal Profile for D'Alessio Disk Temperature",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


def tick_function(X):
    V = 0.1*(X/1000)**(-7/3)
    return V #["%.3f" % z for z in V]



#secax = ax.secondary_xaxis('top', functions = (T_disk, tick_function))
#secax.set_xlabel('Temperature (K)')

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"Isothermal-D'Alessio_GCR for sma and M", dpi=300, bbox_inches='tight')
plt.show()


def rho_nebular(M, a):
    T = T_disk(a)

    P_bot = Ne_reserve * S_Ne / x_Ne / M / M_earth * 1e6
    # P_Ne = 0.12171882539682544 * 1e5
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T) ** 2  # cm

    R_out = [0]
    R_out = min(R_H, R_B)
    R_out_E = R_out / R_E

    # print(R_out)

    A = G * M / Cs_disk(T) ** 2  # CHECK UNITS [cm]

    # P_bot = 10**6#Ne_reserve/S_Ne/x_Ne/M
    R_bot = R_E * M ** 0.25
    exp = np.exp(-A / R_out + A / R_bot)
    rho_bot = P_bot / Cs_disk(T) ** 2  # [g/cm^3]
    rho_ne = rho_bot / exp
    return rho_ne


def sig_ne(a):
    K = np.sqrt(G * M_sun / (a * autocm) ** 3)
    H = Cs_disk(a) / K
    mini_sig = 1.7 * 10 ** 3 * (a * autocm) ** (-3 / 2)
    return H * mini_sig


N = 100
sma_T = sma
M_co = np.logspace(-1, 0, N) * M_Earth
Z_rho = np.zeros((N, N))

# rho_mini = [0]*len(sma_T)
rho_mini = sig_ne(sma_T)
print(rho_mini)

# print(Cs_disk(T_disk(sma_T)))

for i in range(len(M_co)):
    for j in range(len(sma_T)):
        Z_rho[i][j] = rho_nebular(M_co[i], sma_T[j])

# print(Z_rho)

X_rho, Y_rho = np.meshgrid(M_co / M_Earth, sma_T)

fig, ax = plt.subplots(figsize=(10, 8), dpi=80)

origin = 'lower'

min_rho = -8
max_rho = -4

lev_exp = np.arange(min_rho, max_rho, 0.5)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(Y_rho, X_rho, Z_rho, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Semi-Major Axis (AU)')
ax.set_ylabel('Core Mass (M_earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_rho, max_rho, 0.2)  # [-9, -5, -4, -3, -2, -1.5, -1,0.,0.5, 1, 5, 7, 9]
levss = np.power(10., ee)

CS4 = ax.contour(Y_rho, X_rho, Z_rho, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title('Contour Plot of Density to Semi-major Axis and Core Mass')
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)

ax.plot()

#secax = ax.secondary_xaxis('top', functions = (T_disk, tick_function))
#secax.set_xlabel('Temperature (K)')

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Nebular Density (g/cm^3)', rotation=270)
plt.savefig("Density for a and M", dpi=300, bbox_inches='tight')
plt.show()
