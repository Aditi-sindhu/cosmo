
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

c = 3e5  
H0 = 70.0

def Hubble_parameter(z, Omega_m, Omega_r, Omega_DE, w_DE):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_DE * np.exp(3 * quad(lambda z_prime: (1 + w_DE) / (1 + z_prime), 0, z)[0]))

def luminosity_distance(z, Omega_m, Omega_r, Omega_DE, w_DE):
    integrand = lambda z_prime: 1 / Hubble_parameter(z_prime, Omega_m, Omega_r, Omega_DE, w_DE)
    result, _ = quad(integrand, 0, z)
    return (1 + z) * c * result

def comoving_distance(z, Omega_m, Omega_r, Omega_DE, w_DE):
    integrand = lambda z_prime: 1 / Hubble_parameter(z_prime, Omega_m, Omega_r, Omega_DE, w_DE)
    result, _ = quad(integrand, 0, z)
    return c * result

def angular_diameter_distance(z, Omega_m, Omega_r, Omega_DE, w_DE):
    return comoving_distance(z, Omega_m, Omega_r, Omega_DE, w_DE) / (1 + z)

def proper_distance(z, Omega_m, Omega_r, Omega_DE, w_DE):
    a_t0 = 1  
    return a_t0 * comoving_distance(z, Omega_m, Omega_r, Omega_DE, w_DE) / (1+z)

z_values = np.linspace(0, 2, 1000)
#z_values = [0, 2, 6, 1100]

# Flat Lambda Cold Dark Matter
Omega_m1 = 0.3
Omega_r1 = 9.47e-5  
Omega_DE1 = 0.7
w_DE1 = -1
DL1 = [luminosity_distance(z, Omega_m1, Omega_r1, Omega_DE1, w_DE1) for z in z_values]
DC1 = [comoving_distance(z, Omega_m1, Omega_r1, Omega_DE1, w_DE1) for z in z_values]
DA1 = [angular_diameter_distance(z, Omega_m1, Omega_r1, Omega_DE1, w_DE1) for z in z_values]
DP1 = [proper_distance(z, Omega_m1, Omega_r1, Omega_DE1, w_DE1) for z in z_values]

# Flat dark energy (w_DE = -0.7) Cold Dark Matter
Omega_m2 = 0.3
Omega_r2 = 9.47e-5
Omega_DE2 = 0.7
w_DE2 = -0.7
DL2 = [luminosity_distance(z, Omega_m2, Omega_r2, Omega_DE2, w_DE2) for z in z_values]
DC2 = [comoving_distance(z, Omega_m2, Omega_r2, Omega_DE2, w_DE2) for z in z_values]
DA2 = [angular_diameter_distance(z, Omega_m2, Omega_r2, Omega_DE2, w_DE2) for z in z_values]
DP2 = [proper_distance(z, Omega_m2, Omega_r2, Omega_DE2, w_DE2) for z in z_values]

# Flat dark energy (w_DE = -1.2) Cold Dark Matter
Omega_m3 = 0.3
Omega_r3 = 9.47e-5
Omega_DE3 = 0.7
w_DE3 = -1.2
DL3 = [luminosity_distance(z, Omega_m3, Omega_r3, Omega_DE3, w_DE3) for z in z_values]
DC3 = [comoving_distance(z, Omega_m3, Omega_r3, Omega_DE3, w_DE3) for z in z_values]
DA3 = [angular_diameter_distance(z, Omega_m3, Omega_r3, Omega_DE3, w_DE3) for z in z_values]
DP3 = [proper_distance(z, Omega_m3, Omega_r3, Omega_DE3, w_DE3) for z in z_values]

# Flat Lambda Cold Dark Matter with Ωm = 0.8
Omega_m4 = 0.8
Omega_r4 = 9.47e-5
Omega_DE4 = 1 - Omega_m4
w_DE4 = -1
DL4 = [luminosity_distance(z, Omega_m4, Omega_r4, Omega_DE4, w_DE4) for z in z_values]
DC4 = [comoving_distance(z, Omega_m4, Omega_r4, Omega_DE4, w_DE4) for z in z_values]
DA4 = [angular_diameter_distance(z, Omega_m4, Omega_r4, Omega_DE4, w_DE4) for z in z_values]
DP4 = [proper_distance(z, Omega_m4, Omega_r4, Omega_DE4, w_DE4) for z in z_values]

# Flat Lambda Cold Dark Matter with T0 = 100 K
Omega_m5 = 0.3
Omega_r5 = 9.47e-5
Omega_DE5 = 0.7
w_DE5 = -1
DL5 = [luminosity_distance(z, Omega_m5, Omega_r5, Omega_DE5, w_DE5) for z in z_values]
DC5 = [comoving_distance(z, Omega_m5, Omega_r5, Omega_DE5, w_DE5) for z in z_values]
DA5 = [angular_diameter_distance(z, Omega_m5, Omega_r5, Omega_DE5, w_DE5) for z in z_values]
DP5 = [proper_distance(z, Omega_m5, Omega_r5, Omega_DE5, w_DE5) for z in z_values]


plt.plot(z_values, DL1, label='Luminosity_distance Flat Lambda CDM')
plt.plot(z_values, DC1, label='comoving_distance Flat Lambda CDM')
plt.plot(z_values, DA1, label='angular_diameter_distance Flat Lambda CDM')
plt.plot(z_values, DP1, label='proper_distance Flat Lambda CDM')

plt.plot(z_values, DL2, label='Luminosity_distance Flat dark energy(w_DE = -0.7) CDM')
plt.plot(z_values, DC2, label='comoving_distance Flat dark energy(w_DE = -0.7) CDM')
plt.plot(z_values, DA2, label='angular_diameter_distance Flat dark energy(w_DE = -0.7) CDM')
plt.plot(z_values, DP2, label='proper_distance Flat dark energy(w_DE = -0.7) CDM')

plt.plot(z_values, DL3, label='Luminosity_distance Flat dark energy (w_DE = -1.2) CDM')
plt.plot(z_values, DC3, label='comoving_distance Flat dark energy (w_DE = -1.2) CDM')
plt.plot(z_values, DA3, label='angular_diameter_distance Flat dark energy (w_DE = -1.2) CDM')
plt.plot(z_values, DP3, label='proper_distance Flat dark energy (w_DE = -1.2) CDM')

plt.plot(z_values, DL4, label='Luminosity_distance Flat Lambda CDM with Ωm = 0.8')
plt.plot(z_values, DC4, label='comoving_distance Flat Lambda CDM with Ωm = 0.8')
plt.plot(z_values, DA4, label='angular_diameter_distance Flat Lambda CDM with Ωm = 0.8')
plt.plot(z_values, DP4, label='proper_distance Flat Lambda CDM with Ωm = 0.8')

plt.plot(z_values, DL5, label='Luminosity_distance Flat Lambda CDM with T0 = 100 K')
plt.plot(z_values, DC5, label='comoving_distance Flat Lambda CDM with T0 = 100 K')
plt.plot(z_values, DA5, label='angular_diameter_distance Flat Lambda CDM with T0 = 100 K')
plt.plot(z_values, DP5, label='proper_distance Flat Lambda CDM with T0 = 100 K')

#plt.legend()
plt.show()
