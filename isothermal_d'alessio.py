from scipy import integrate
import scipy.integrate as integrate
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker, cm, colors
from scipy.interpolate import interp1d




SMALL_SIZE = 13
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}

M_Earth = 5.972e+27
mu = 2.34
m_p = 1.67 * 10 ** (-24)  # [g]
k_B = 1.38 * 10 ** (-16)  # [erg/K]
G = 6.67408 * 10 ** (-8)  # [dyne cm^2/g^2]
M_mars = 0.64171 * 10 ** (27)  # [g]
M_sun = 1.989 * 10 ** (33)  # [g]
M_earth = 5.98 * 10 ** 27  # [g]
R_E = 6.378 * 10 ** 8  # [cm]
autocm = 1.496 * 10 ** (13)
Ne_reserve = 64.5e15
S_Ne = 2.7e-12
x_Ne = 2.1e-5

T_0 = 1000


# define temperature
def T_disk(dis):
    a = np.loadtxt(open("Dalessio-1.csv", "rb"), delimiter=",", skiprows=1)
    b = np.loadtxt(open("Dalessio-2.csv", "rb"), delimiter=",", skiprows=1)

    R_1 = np.sort(a[:, 0])
    T_1 = a[:, 1]

    R_2 = np.sort(b[:, 0])
    T_2 = b[:, 1]

    f1 = interp1d(R_1, T_1, kind='cubic', fill_value='extrapolate')
    f2 = interp1d(R_2, T_2, kind='cubic', fill_value='extrapolate')
    x_1 = np.linspace(np.min(R_1), 0.1, 200, endpoint=True)
    x_2 = np.linspace(0.1, np.max(R_2), 200, endpoint=True)

    if dis < 0.1:
        return f1(dis)

    else:
        return f2(dis)

    #return T_0 * (a / 0.1) ** (-3 / 7)  # [K]


def Cs_disk(T):
    return np.sqrt(k_B * T / mu / m_p)  # CHECK UNITS [cm/s]


def M_int(x, a, M, P_bot, T):
    # print(Cs_disk(T))
    # compute the R_out based on Boundy/Hill Radisu
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T) ** 2  # cm

    R_out = [0]
    R_out = min(R_H, R_B)
    R_out_E = R_out / R_E

    A = G * M / Cs_disk(T) ** 2  # CHECK UNITS [cm]

    # P_bot = 10**6#Ne_reserve/S_Ne/x_Ne/M
    R_bot = R_E * M ** 0.25
    exp = np.exp(-A / R_out + A / R_bot)
    rho_bot = P_bot / Cs_disk(T) ** 2  # [g/cm^3]
    rho_ne = rho_bot / exp

    return 4 * np.pi * rho_ne * (x ** 2 * np.exp(A / R_bot - A / x))


M_array = np.linspace(1e-1, 10, 5000)
M_P_Constant = [0] * len(M_array)
M_P_Ne_limit = [0] * len(M_array)
M_P_10bar = [0] * len(M_array)
M_P_5bar = [0] * len(M_array)


def r_max(M, T):
    return G * M * M_Earth / Cs_disk(T) ** 2

T_surf = 2000 #[K] using a quasi-isothermal profile

def GCR_isothermal(a, M):
    T = T_surf
    r_min = R_E * (M) ** 0.25
    # P_Ne = Ne_reserve/S_Ne/x_Ne/M/M_Earth
    P_Ne = 0.12171882539682544 * 1e5
    res = quad(M_int, r_min, r_max(M, T), args=(a, M * M_Earth, P_Ne, T))[0]
    return res / M / M_Earth


N = 100
sma = np.logspace(-1, 1, N)
M_co = np.logspace(-1, 0, N)
Z = np.zeros((N, N))

X, Y = np.meshgrid(M_co, sma)
for i in range(len(M_co)):
    for j in range(len(sma)):
        Z[i][j] = GCR_isothermal(M_co[i], sma[j])

fig , ax = plt.subplots(figsize=(10, 8), dpi=80)


max_GCR = 2
min_GCR = -8

origin = 'lower'
lev_exp = np.arange(min_GCR, max_GCR, 0.2)
levs = np.power(10, lev_exp)
# print(lev_exp)
cs = ax.contourf(Y, X, Z, levs, norm=colors.LogNorm(), cmap='Spectral')
# cs = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap=cm.PuBu_r)
ax.set_xlabel('Semi-major axis(AU)')
ax.set_ylabel('Core Mass (M_Earth)')
ax.set_yscale('log')
ax.set_xscale('log')
ee = np.arange(min_GCR, max_GCR, 1)
levss = np.power(10., ee)

CS4 = ax.contour(Y, X, Z, levss,
                 colors=('cyan',),
                 linewidths=(1,),
                 origin=origin)
ax.set_title(f'Quasi-Isothermal Profile for DAlesson Disk Temperature, T = {T_surf}',fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


def tick_function(X):
    V = 0.1*(X/1000)**(-7/3)
    return V #["%.3f" % z for z in V]



#secax = ax.secondary_xaxis('top', functions = (T_disk, tick_function))
#secax.set_xlabel('Temperature (K)')

cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel('Gas to Core Ratio (GCR)', rotation=270)
plt.savefig(f"quasi-D'Alesson_GCR for sma and M, T = {T_surf}", dpi=300, bbox_inches='tight')
plt.show()


def rho_nebular(M, a):
    T = T_disk(a)

    P_bot = 0.0082171882539682544 * 1e5
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
ax.set_ylabel('Core Mass (M_Earth)')
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
