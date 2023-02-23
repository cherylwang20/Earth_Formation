import matplotlib.pyplot as plt
import numpy as np

from functions_const import *

print(M_earth/M_sun)
#we first write down the orbital crossing timescale under

print(Cs_disk(T_disk_chiang(0.1))/np.sqrt(G*0.5*M_earth/0.1/autocm))

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
    return rho*6*10**(-6)/MMSN_V(a)

a_semi = np.linspace(0.1, 1, 100)
M_embryos = 0.5

ecc = [0, 1E-3, 1E-2]

fig , ax = plt.subplots(figsize=(9, 7), dpi=80)

for i in range(len(ecc)):
    plt.plot(a_semi, max_rho(a_semi, M_embryos, ecc[i], tx(M_embryos,M_sun,a_semi,ecc[i])), linewidth =2, label = f'{ecc[i]}')
plt.yscale('log')
#plt.xscale('log')
applyplotstyle()
plt.ylabel(r'$\rho_{neb}$/$\rho_{MMSN}$')
plt.xlabel('Semimajor-Axis (AU)')
plt.title(f'M = {M_embryos}$M_\oplus$')
plt.xlim([0.1,1])
plt.ylim([10**-8, 10**-6])
plt.legend(title = 'eccentricity (e)', bbox_to_anchor=(0.9, 0.905))
plt.gca().invert_yaxis()
plt.savefig(f"Nebular_density_tx_tdamp_2", dpi=300, bbox_inches='tight')
plt.show()
plt.close()

