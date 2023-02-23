import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

print(M_earth/M_sun)
#we first write down the orbital crossing timescale under

def tx(M_p, M_star,a, e):
    M_p = 2*M_p
    R_mH = (2*M_p/3/M_star)**(1/3)*a*autocm #accounting for the masses of two embryos
    print(R_mH)
    del_a = 20*R_mH #change accordingly
    k = del_a/R_mH
    h = k/2*(2*M_p/3/M_star)**(1/3)
    A = -2 + e/h - 0.27* np.log10(M_p/M_star)
    B = 18.7  - 16.8*e/h + (1.1 - 1.2*e/h)*np.log10(M_p/M_star)
    return 10**(A + B*np.log10(k/2.3))*a**1.5 #in years

#t_damp = 0.1T_x
def max_rho(a,M_p,t_x):
    v_rel = Cs_disk(T_disk_chiang(a))
    t_x = t_x*3.154e+7
    return v_rel**3/(4*np.pi*G**2*M_p*0.1*t_x)

#let the embryos be 0.5 earth masses
# assuming eccentricity for 0.01
M_embryos = 0.5*M_earth
eccen = 0.115/2
semi = 0.1 #AU

ec = np.logspace(-4, -0.5, 100)
M_e = np.linspace(0.1, 1, 100)*M_earth
a_semi = np.linspace(0.1, 1, 100)

t_embryos = tx(M_embryos, M_sun, semi, eccen)
rho_embryos = max_rho(semi, M_embryos, t_embryos)

print(rho_embryos)

plt.plot(ec, max_rho(semi, M_embryos, tx(M_embryos,M_sun,semi,ec))/MMSN_V(semi))
plt.yscale('log')
plt.xscale('log')
plt.ylim([10**-10,1])
plt.gca().invert_yaxis()
plt.show()
plt.close()

fig , ax = plt.subplots(figsize=(9, 7), dpi=80)
plt.plot(a_semi, max_rho(a_semi, M_embryos, tx(M_embryos,M_sun,a_semi,0.1)), linewidth =2, label = f'0.1')
plt.plot(a_semi, max_rho(a_semi, M_embryos, tx(M_embryos,M_sun,a_semi,0.05)), linewidth =2, label = f'0.05')
plt.plot(a_semi, max_rho(a_semi, M_embryos, tx(M_embryos,M_sun,a_semi,0.01)), linewidth =2, label = f'0.01')
plt.yscale('log')
plt.xscale('log')
applyplotstyle()
plt.ylabel(f'Nebular Density (g/cm$^3$)')
plt.xlabel('Semimajor-Axis (AU)')
plt.title(f'M = {M_embryos/M_earth}$M_\oplus$')
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.savefig(f"Nebular_density_tx_tdamp", dpi=300, bbox_inches='tight')
plt.show()
plt.close()
