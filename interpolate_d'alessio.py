import numpy as np
from scipy import interpolate
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

a = np.loadtxt(open("Dalessio-1.csv", "rb"), delimiter=",", skiprows=1)
b = np.loadtxt(open("Dalessio-2.csv", "rb"), delimiter=",", skiprows=1)

R_1 = np.sort(a[:,0])
T_1 = a[:,1]

R_2 = np.sort(b[:,0])
T_2 = b[:,1]

f1 = interp1d(R_1, T_1, kind='cubic', fill_value='extrapolate')
f2 = interp1d(R_2, T_2, kind='cubic', fill_value='extrapolate')
x_1 = np.linspace(np.min(R_1), 0.1,200,endpoint=True)
x_2 = np.linspace(0.1, np.max(R_2),200,endpoint=True)

plt.loglog(x_1,f1(x_1))
plt.loglog(x_2,f2(x_2))
plt.show()
