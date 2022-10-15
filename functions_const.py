#for atmospheric accretion
import numpy as np
from scipy import integrate
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib import ticker as mticker
from scipy import interpolate
from matplotlib import ticker, cm, colors

SMALL_SIZE = 15
matplotlib.rc('font', size=SMALL_SIZE, family='Arial')
matplotlib.rc('axes', titlesize=SMALL_SIZE)

hfont = {'fontname': 'Helvetica'}

def applyplotstyle(ax):
    origin = 'lower'
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.xlim([1E-8, 1])
    ax.set_xscale('log')
    plt.ylabel(r'Core Mass (M_$\oplus$)')
    ax.set_xlabel(r'$\rho_{neb}$/$\rho_{MMSN}$')
    ax.xaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax.xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
    plt.tick_params(direction='in', length=6, width=1.5)
    plt.tick_params(axis="y", which='major', direction="in", width=1.5, length=8, pad=4)
    plt.tick_params(axis="y", which='minor', direction="in", width=1.5, length=4, pad=4)
    plt.tick_params(axis="x", which='major', direction="in", width=1.5, length=8, pad=7)
    plt.tick_params(axis="x", which='minor', direction="in", width=1.5, length=4, pad=7)
    plt.tick_params(right=True, top=True, labelright=False, labeltop=False)

mu = 2.34
m_p = 1.67 * 10 ** (-24)  # [g]
k_B = 1.38 * 10 ** (-16)  # [erg/K]
G = 6.67408 * 10 ** (-8)  # [dyne cm^2/g^2]
M_mars = 0.64171 * 10 ** (27)  # [g]
M_sun = 1.989 * 10 ** (33)  # [g]
M_earth = 5.98 * 10 ** 27  # [g]
R_E = 6.378 * 10 ** 8  # [cm]
autocm = 1.496 * 10 ** (13)
Ne_reserve = 2.735 * 10 ** (15)
S_Ne = 2.7e-12
x_Ne = 2.1e-5
sigma = 5.6704 * 10 ** (5) #cgs
m_Ne = 3.35 * 10**(-23)#change to H for testing #g/mol
m_He = 6.64*10**(-24)
m_H = 1.67*10**(-24)
mole_H= 1.00784
mole_He = 4.002602
mole_Ne = 20.1797
sigma_neon = 0.24e-14 #cm^2


# define temperature, which we only use as the disk temperature from Chiang
def T_disk_chiang(a):
    return 150 * (a ) ** (-3 / 7)  # [K]

#define disk temperature from d'alessio, through interporlation
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


#define sounds speed
def Cs_disk(T):
    return np.sqrt(k_B * T / mu / m_p)  # CHECK UNITS [cm/s]

#isothermal atmospheric mass
def M_int(x, a, M, P_bot, T):
    # compute the R_out based on Boundy/Hill Radisu
    R_H = a * (M  / 3 / M_sun) ** (1 / 3) * autocm # cm
    R_B = G * M  / Cs_disk(T) ** 2   # cm


    R_out = min(R_H, R_B)
    R_out_E = R_out / R_E

    A = G * M / Cs_disk(T) ** 2  # CHECK UNITS [cm]

    # P_bot = 10**6#Ne_reserve/S_Ne/x_Ne/M
    R_bot = R_E * M ** 0.25
    exp = np.exp(-A / R_out + A / R_bot)
    rho_bot = P_bot / Cs_disk(T) ** 2  # [g/cm^3], using the qusi-isothermal approximation.
    rho_ne = rho_bot / exp

    return 4 * np.pi * rho_ne * (x ** 2 * np.exp(A / R_bot - A / x))


def quasi_M(x, a, M, P_bot, T_bot, T_top):
    # compute the R_out based on Boundy/Hill Radisu
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm

    R_out = min(R_H, R_B)
    R_out_E = R_out / R_E

    A = G * M / Cs_disk(T_bot) ** 2  # CHECK UNITS [cm]

    # P_bot = 10**6#Ne_reserve/S_Ne/x_Ne/M
    R_bot = R_E * M ** 0.25
    exp = np.exp(-A / R_out + A / R_bot)
    rho_bot = P_bot / Cs_disk(T_bot) ** 2  # [g/cm^3], using the qusi-isothermal approximation.
    rho_ne = rho_bot / exp

    return 4 * np.pi * rho_ne * (x ** 2 * np.exp(A / R_bot - A / x))


def r_max(M, T):
    return G * M * M_earth / Cs_disk(T) ** 2

def r_min(M):
    return R_E * (M) ** 0.25

#GCR for an isothermal profile
def GCR_isothermal(a, M, T_fun):
    T = T_fun(a)
    r_min = R_E * (M) ** 0.25
    S_Ne = np.exp(-lnk_Ne_Earth(T, 0.1)[0])
    P_Ne = Ne_reserve * S_Ne / x_Ne / M / M_earth * 1e6
    res = quad(M_int, r_min, r_max(M, T), args=(a, M * M_earth, P_Ne, T))[0] #Mass given in grams
    return res / M / M_earth


def GCR_quasi(a, M, T_fun, T_bot):
    T_top = T_fun(a)
    print(M)
    r_min = R_E * (M) ** 0.25
    S_Ne = np.exp(-lnk_Ne_Earth(T_bot, 0.1)[0])
    print(S_Ne)
    P_Ne = Ne_reserve * S_Ne / x_Ne / M / M_earth * 1e6
    print(P_Ne)
    res = quad(quasi_M, r_min, r_max(M, T_bot), args=(a, M * M_earth, P_Ne, T_bot, T_top))[0] #Mass given in grams
    return res / M / M_earth

#Henry's law

alpha_Ne = 0.401
beta_Ne = -43.36
nu_ca = 13.30
nu_ca_MgO = 7.57
nu_s = -7.18
nu_s_MgO = -2.2
lambda_sio2 = -770.13
lambda_MgO = 2118.9
k_sio2 = -2.06e-6
k_MgO = 4.7e-5
V_Ne = 26.5
Mol_SiO2 = 60.09
Mol_MgO = 40.3
T_ref = 1300  # degree
P_ref = 0.1  # MPa
R_gas = 8.314 * 10**7 #erg/mol/K
av = 6.02214076 * 10**23 #mol^-1


def lnk_Ne_Earth(T, P):
    mol_sio2 = 0.446
    mol_MgO = 0.554
    fun = nu_ca*mol_sio2 + nu_ca_MgO*mol_MgO + nu_s + (lambda_sio2* mol_sio2 + mol_MgO*lambda_MgO) * (1 / T - 1 / T_ref) + (k_sio2*mol_sio2 + k_MgO*mol_MgO) * (P - P_ref) * 10
    result = alpha_Ne * (100 - 100 / (49.14 * 112.33 / 271.85) * fun) + beta_Ne
    return result, fun

#accretion luminosity
def lum_fun(M, a, M_rate):
    Rp = R_E * M ** 0.25
    return G*M*M_earth*M_rate(M,a)/Rp

def M_adi(x, gamma, M, a, T_top, rho_neb):
    B = (gamma - 1) / gamma * G * M / Cs_disk(T_top) ** 2
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm
    R_out = min(R_H, R_B)
    r_min = R_E * (M / M_earth) ** 0.25
    fun = rho_neb*(1 + B*(1/x - 1/R_out))**(1/(gamma - 1))
    return 4 * np.pi * x**2 * fun

def adi_P_T(gamma, M, a, T_top, rho_neb):
    B = (gamma - 1) / gamma * G * M / Cs_disk(T_top) ** 2
    R_H = a * (M / 3 / M_sun) ** (1 / 3) * autocm  # cm
    R_B = G * M / Cs_disk(T_top) ** 2  # cm
    R_out = min(R_H, R_B)
    #print(R_out)
    r_min = R_E * (M / M_earth) ** 0.25
    rho_bot = rho_neb*(1 + B*(1/r_min - 1/R_out))**(1/(gamma - 1))
    C = (gamma - 1) / gamma * mu / R_gas * G * M
    T_bot = T_top - C/R_out + C/r_min
    P_bot = rho_bot*Cs_disk(T_bot)
    return rho_bot, P_bot, R_out, T_bot

def MMSN_V(a):
    MMSN_disk = 1700* (a)**(-3./2.)

    K = np.sqrt(G*M_sun/(a*autocm)**3) #[1/s]
    H = Cs_disk(T_disk_chiang(a))/K

    return MMSN_disk/H

def tx(M_p, M_star,a, e):
    M_p = 2*M_p*M_earth
    R_mH = (2*M_p/3/M_star)**(1/3)*a*autocm #accounting for the masses of two embryos
    #print(R_mH)
    del_a = 16*R_mH #change accordingly
    k = del_a/R_mH
    h = k/2*(2*M_p/3/M_star)**(1/3)
    A = -2 + e/h - 0.27* np.log10(M_p/M_star)
    B = 18.7  - 16.8*e/h + (1.1 - 1.2*e/h)*np.log10(M_p/M_star)
    return 10**(A + B*np.log10(k/2.3))*a**1.5 #in years

#t_damp = 0.1T_x
def max_rho(a,M_p,e,t_x):
    v_ratio = Cs_disk(T_disk_chiang(a))/np.sqrt(G*M_p*M_earth/a/autocm)
    rho = 0.5*(T_disk_chiang(a)/1000)**1.5/M_p*(1 + 0.25*(e/v_ratio)**3)/0.1/t_x
    return rho*6*10**(-6)

def des_rho(total_neon, degas):
    N = 100
    M_co = np.linspace(0.1, 1, N)
    rho = np.logspace(-18, -5, N)

    f = interpolate.interp2d(M_co, rho, total_neon, kind='cubic')
    f3 = interpolate.interp1d(M_co, degas)

    des_rho = [0] * len(M_co)
    rho_den = np.logspace(-18, -5, 10 ** 5)

    for j in range(len(M_co)):
        for i in range(10 ** 5):
            if np.abs((f(M_co[j], rho_den[i]) - f3(M_co[j]))) / f3(M_co[j]) < 0.001:
                print(f3(M_co[j]), f(M_co[j], rho_den[i]))
                print(rho_den[i])
                des_rho[j] = rho_den[i]
                break
    return des_rho

