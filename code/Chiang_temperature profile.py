from functions_const import *

def T_disk(a):
    return 150 * (a ) ** (-3 / 7)

a = np.loadtxt(open("Dalessio-1.csv", "rb"), delimiter=",", skiprows=1)
b = np.loadtxt(open("Dalessio-2.csv", "rb"), delimiter=",", skiprows=1)

R_1 = np.sort(a[:,0])
T_1 = a[:,1]

R_2 = np.sort(b[:,0])
T_2 = b[:,1]

f1 = interp1d(R_1, T_1, kind='cubic', fill_value='extrapolate')
f2 = interp1d(R_2, T_2, kind='cubic', fill_value='interpolation')
#x_1 = np.linspace(0.4, 0.1,200,endpoint=True)
x_2 = np.linspace(0.1, 10,100,endpoint=True)
fig, ax = plt.subplots(figsize=(10, 8), dpi=80)
#plt.loglog(x_1,f1(x_1),  label = 'D-Alessio', c = 'blue')
plt.loglog(x_2,f2(x_2), label = 'D-Alessio', c = 'blue')


a = np.linspace(0.1, 10, 500)
plt.loglog(a, T_disk(a), label = 'Chiang', c = 'red')

plt.legend()
plt.xlabel('Semi-Major Axis (AU)')
plt.ylabel('Temperature (K)')
plt.savefig(f"Chiang vs D'Alessio", dpi=300, bbox_inches='tight')

plt.show()

M = np.linspace(0.1, 1, 5)

T = T_disk(a)
T = f2(a)

def Cs_disk(T):
    return np.sqrt(k_B * T / mu / m_p)  # CHECK UNITS [cm/s]

def radius(M):
    R_H = a * (M * M_earth/ 3 / M_sun) ** (1 / 3) * autocm / R_E  # cm
    R_B = G * M * M_earth / Cs_disk(T) ** 2 / R_E  # cm
    return R_H, R_B

fig, ax = plt.subplots(figsize=(10, 8), dpi=80)

for i in range(len(M)):
    plt.loglog(a, radius(M[i])[0], label=('' if i==0 else '_')  + 'Hill Radius', c = 'orange')
    plt.loglog(a, radius(M[i])[1], label=('' if i==0 else '_')  +'Bondi Radius', c = 'green')
plt.xlabel('Semi-Major Axis (AU)')
plt.ylim([1, 3000])
plt.ylabel('Radius (R_E)')
plt.legend()

plt.savefig(f"D'alessio-Hill vs Bondi Radius", dpi=300, bbox_inches='tight')
plt.show()