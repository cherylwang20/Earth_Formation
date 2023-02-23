from functions_const import *

M_t = 0.5
P = 10E-1
T = 2500
T_now = 1500

def Ne_res(M_core, P, T):
    # Ne_con = 3.42e-14
    S_Ne = np.exp(-lnk_Ne_Earth(T, P)[0])
    # print(S_Ne)
    return M_core * P / S_Ne * x_Ne * 5.972e+24 * 1000  # g

N1 = Ne_res(M_t, P, T)
N2 = Ne_res(M_t, P, T_now)

print((N1 - N2)/N1)