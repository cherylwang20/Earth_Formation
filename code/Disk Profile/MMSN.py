import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from functions_const import *

m_v = 0.857*M_earth
rmhh = ((M_earth + m_v)/3/M_sun)**(1/3)*0.86

kk = (1-0.723)/rmhh

Rmh = (2*M_earth/3/M_sun)**(1/3)*0.86
k = 0.28/Rmh

T_vir = G*M_earth*2.37*m_p/3/r_min(1)/k_B
print(T_vir)

mu_g = 8*10**6*k_B*2.5*R_sun/(G*0.5*M_sun)
print(mu_g/m_p)

N = 100
a = np.linspace(0.1, 1, N)
MMSN_disk = 1700* (a)**(-3./2.)

K = np.sqrt(G*M_sun/(a*autocm)**3) #[1/s]
H = Cs_disk(T_disk_chiang(a))/K

MMSN_v = MMSN_disk/H

print(1700/autocm**(3/2))

#plt.plot(a, MMSN_v)
plt.plot(a, 1700*np.sqrt(G*M_sun)/(Cs_disk(T_disk_chiang(a))*a**3*autocm**(3/2)))
plt.yscale('log')
#applyplotstyle()
plt.show()
plt.close()

M_co = np.linspace(0.1, 1, N)
rho = np.logspace(-18, -5, N)

total_neon = np.loadtxt(open("output_9.csv"), delimiter=",")
degas_99 = np.loadtxt(open("degas_9.csv"),delimiter=",")

print(degas_99)
#print(total_neon)
#print(degas_99)

f = interpolate.interp2d(M_co, rho, total_neon, kind='cubic')
f3 = interpolate.interp1d(M_co,degas_99)

des_rho = [0]*len(M_co)
rho_den = np.logspace(-18,-5, 10**5)

for j in range(len(M_co)):
    for i in range(10**5):
        if np.abs((f(M_co[j], rho_den[i])-f3(M_co[j])))/f3(M_co[j]) < 0.001:
            print(f3(M_co[j]),f(M_co[j], rho_den[i]))
            print(rho_den[i])
            des_rho[j] = rho_den[i]
            break

print(f'the density is: {des_rho}')
np.savetxt('desire_rho_9.csv',des_rho,delimiter=",")