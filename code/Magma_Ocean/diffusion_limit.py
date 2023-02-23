import matplotlib.pyplot as plt

from functions_const import *

N = 100
M_co = np.linspace(0.1, 1, N)
T = 10000;

r_sonic = 0.5*r_max(M_co, T)
sigma_Ne_H = 1.34*10**(-15)

n_photo = 10**5
scale_height = k_B*T*r_sonic**2/G/M_co/M_earth/m_Ne
mfp = 1/sigma_Ne_H/n_photo #photoionization n = 10**9; for sonic point: 10**7

print(scale_height)

ratio = mfp/scale_height

fig, ax = plt.subplots(figsize=(9,7), dpi=80,constrained_layout=True)
plt.plot(M_co, ratio)
plt.xscale('log')
plt.xlabel('Core Mass (M_E)')
plt.ylabel('Ratio of mfp to Scale Height')
plt.title(f'n = {n_photo:.2e}')
#plt.yscale('log')
plt.show()

