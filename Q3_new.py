
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math
H0 = 0.07    #in Gyr inv

def age_of_universe(z, Omega_m, Omega_r, Omega_DE, w_DE):
    result, _ = quad(lambda z_prime: age_integrand(z_prime, Omega_m, Omega_r, Omega_DE, w_DE), z, 8000)
    return result

def age_integrand(z, Omega_m, Omega_r, Omega_DE, w_DE):
    return 1 / ((1 + z) * Hubble_parameter(z, Omega_m, Omega_r, Omega_DE, w_DE))

def Hubble_parameter(z, Omega_m, Omega_r, Omega_DE, w_DE):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_r * (1 + z)**4 + Omega_DE * (1+z)**(-3*(1+w_DE)))


z_values = [0, 2, 6, 1100]

#plt.figure(figsize=(12, 8))

# Flat Lambda Cold Dark Matter
Omega_m1 = 0.3
Omega_r1 = 9.47e-5  
Omega_DE1 = 0.7
w_DE1 = -1
age1 = [age_of_universe(z, Omega_m1, Omega_r1, Omega_DE1, w_DE1) for z in z_values]

# Flat dark energy (w_DE = -0.7) Cold Dark Matter
Omega_m2 = 0.3
Omega_r2 = 9.47e-5
Omega_DE2 = 0.7
w_DE2 = -0.7
age2 = [age_of_universe(z, Omega_m2, Omega_r2, Omega_DE2, w_DE2) for z in z_values]

# Flat dark energy (w_DE = -1.2) Cold Dark Matter
Omega_m3 = 0.3
Omega_r3 = 9.47e-5
Omega_DE3 = 0.7
w_DE3 = -1.2
age3 = [age_of_universe(z, Omega_m3, Omega_r3, Omega_DE3, w_DE3) for z in z_values]

# Flat Lambda Cold Dark Matter with Ωm = 0.8
Omega_m4 = 0.8
Omega_r4 = 9.47e-5
Omega_DE4 = 1 - Omega_m4
w_DE4 = -1
age4 = [age_of_universe(z, Omega_m4, Omega_r4, Omega_DE4, w_DE4) for z in z_values]

# Flat Lambda Cold Dark Matter with T0 = 100 K
Omega_m5 = 0.3
Omega_r5 = 9.47e-5
Omega_DE5 = 0.7
w_DE5 = -1
age5 = [age_of_universe(z, Omega_m5, Omega_r5, Omega_DE5, w_DE5) for z in z_values]

print(age1)
print(age2)
print(age3)
print(age4)
print(age5)
plt.plot(z_values, age1, label='Flat Lambda CDM')
plt.plot(z_values, age2, label='Flat DE (w_DE = -0.7) CDM')
plt.plot(z_values, age3, label='Flat DE (w_DE = -1.2) CDM')
plt.plot(z_values, age4, label='Flat Lambda CDM (Ωm = 0.8)')
plt.plot(z_values, age5, label='Flat Lambda CDM (T0 = 100 K)')

plt.title('Age of the Universe at Different Redshifts')
plt.xlabel('Redshift')
plt.ylabel('Age of the Universe (Gyr)')
plt.legend()
plt.show()
