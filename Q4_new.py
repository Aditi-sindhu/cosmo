
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

H0 = 70.0  

def hubble_integrand(z, Omega_m, Omega_DE, w_DE):
    return 1 / np.sqrt(Omega_m * (1 + z)**3 + Omega_DE * np.exp(3 * quad(lambda z_prime: (1 + w_DE) / (1 + z_prime), 0, z)[0]))

def hubble_parameter(z, Omega_m, Omega_DE, w_DE):
    return H0 * np.sqrt(Omega_m * (1 + z)**3 + Omega_DE * np.exp(3 * quad(lambda z_prime: (1 + w_DE) / (1 + z_prime), 0, z)[0]))

def lookback_time_integrand(z, Omega_m, Omega_DE, w_DE):
    return 1 / ((1 + z) * hubble_parameter(z, Omega_m, Omega_DE, w_DE))

def lookback_time(z, Omega_m, Omega_DE, w_DE):
    result, _ = quad(lambda z_prime: lookback_time_integrand(z_prime, Omega_m, Omega_DE, w_DE), 0, z)
    return result

z_values = [0, 2, 6, 1100]

plt.figure(figsize=(12, 8))

# Flat Lambda Cold Dark Matter
Omega_m1 = 0.3
Omega_DE1 = 0.7
w_DE1 = -1
t1 = [lookback_time(z, Omega_m1, Omega_DE1, w_DE1) for z in z_values]

# Flat dark energy (w_DE = -0.7) Cold Dark Matter
Omega_m2 = 0.3
Omega_DE2 = 0.7
w_DE2 = -0.7
t2 = [lookback_time(z, Omega_m2, Omega_DE2, w_DE2) for z in z_values]

# Flat dark energy (w_DE = -1.2) Cold Dark Matter
Omega_m3 = 0.3
Omega_DE3 = 0.7
w_DE3 = -1.2
t3 = [lookback_time(z, Omega_m3, Omega_DE3, w_DE3) for z in z_values]

#  Flat Lambda Cold Dark Matter with Ωm = 0.8
Omega_m4 = 0.8
Omega_DE4 = 1 - Omega_m4
w_DE4 = -1
t4 = [lookback_time(z, Omega_m4, Omega_DE4, w_DE4) for z in z_values]

# Flat Lambda Cold Dark Matter with T0 = 100 K
Omega_m5 = 0.3
Omega_DE5 = 0.7
w_DE5 = -1
t5 = [lookback_time(z, Omega_m5, Omega_DE5, w_DE5) for z in z_values]

print(t1)
print(t2)
print(t3)
print(t4)
print(t5)
plt.plot(z_values, t1, label='Flat Lambda CDM')
plt.plot(z_values, t2, label='Flat DE (w_DE = -0.7) CDM')
plt.plot(z_values, t3, label='Flat DE (w_DE = -1.2) CDM')
plt.plot(z_values, t4, label='Flat Lambda CDM (Ωm = 0.8)')
plt.plot(z_values, t5, label='Flat Lambda CDM (T0 = 100 K)')

plt.title('Lookback Time of the Universe at Different Redshifts')
plt.xlabel('Redshift')
plt.ylabel('Lookback Time ')
plt.legend()
plt.show()