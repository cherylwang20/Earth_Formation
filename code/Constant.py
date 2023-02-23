class Constant:
    def __init__(self, name, value, units):
        self.name = name
        self.value = value
        self.units = units


mu = Constant('mean molecular weight', 2.37, '')
m_p = Constant('proton mass', 1.67 * 10 ** (-24), 'g')  # [g]
k_B = Constant('Boltzmann Constant', 1.38 * 10 ** (-16) , 'erg/K') # [erg/K]
G = Constant('Gravitatonal Constant',6.67408 * 10 ** (-8),  'dyne cm^2/g^2')
M_mars = 0.64171 * 10 ** (27)  # [g]
M_sun = 1.989 * 10 ** (33)  # [g]
M_earth = 5.98 * 10 ** 27  # [g]
R_E = 6.378 * 10 ** 8  # [cm]
R_sun = 6.957* 10**10 #[cm]
autocm = 1.496 * 10 ** (13)
Ne_reserve = 2.735 * 10 ** (15)
S_Ne = 2.7e-12
x_Ne = 2.1e-5 #molar ratio of Ne in solar composition
sigma = 5.6704 * 10 ** (5) #cgs
m_Ne = 3.35 * 10**(-23)#change to H for testing #g/mol
m_He = 6.64*10**(-24)
m_H = 1.67*10**(-24) #mass of hydrogen
mole_H= 1.00784
mole_He = 4.002602
mole_Ne = 20.1797
sigma_neon = 0.24e-14 #cm^2
tau_disk = 2.5 #Myr
