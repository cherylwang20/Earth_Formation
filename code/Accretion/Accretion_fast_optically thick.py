import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import accretion_fast as a1 #import T_disk, Cs_disk, lum_fast
#import Accretion_Rate as a2
from matplotlib import ticker, cm, colors

SMALL_SIZE = 13
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}

slow_t = 3 * 10**10 #year
fast_t = 3*10**5
yeartos = 3.154e+7
M_earth = 5.98 * 10 ** 27  # [g]
R_E = 6.378 * 10 ** 8  # [cm]
G = 6.67408 * 10 ** (-8)  # [dyne cm^2/g^2]
sigma = 5.6704 * 10 ** (5) #cgs
mu = 2.34
m_p = 1.67 * 10 ** (-24)  # [g]
k_B = 1.38 * 10 ** (-16)  # [erg/K]

xi = 1
alpha = 0
beta = 0
Del_0 = (1 + alpha)/(4 + beta)
k_0 = 0.1 # cm^2/g



def Ra(M, a):
    T_0 = a1.T_disk(a)
    R_B = G * M * M_earth / a1.Cs_disk(T_0) ** 2
    lambda_thick = 1.7 * 10 ** 9 * a ** (11 / 4) / (k_0 / 0.1)
    R_L = (a1.lum_fast(M, a)/(15*np.pi*sigma*T_0**4))**(1/2)
    return R_B*(1+alpha)/np.log(R_B*lambda_thick/R_L**2), R_B

N = 100
sma = np.logspace(-1, 1, N)
Mp = np.logspace(-1, 0, N)

plt.plot(Mp, Ra(Mp,1)[0])
plt.plot(Mp, np.zeros(100)*R_E*(Mp**0.25))
plt.show()


def T(M,a):
    T_0 = a1.T_disk(a)
    R_a, R_B = Ra(M, a)
    Rp = R_E * M ** 0.25
    return T_0*(xi + Del_0*(R_B/Rp - R_B/R_a))



plt.plot(Mp, T(Mp, 1))
plt.show()

X, Y = np.meshgrid(Mp, sma)
Z = T(X, Y)

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
ax.set_title(f"T_surf Optically Thick in Fast regime",fontsize=16)
ax.clabel(CS4, fmt='%1.1e', colors='black', fontsize=12)


cbar = fig.colorbar(cs)
cbar.ax.get_yaxis().labelpad = 25
cbar.ax.set_ylabel(r'T_surf', rotation=270)
plt.savefig(f"T_surf_thick_fast", dpi=300, bbox_inches='tight')
plt.show()