
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import h, c, k, pi

T0 = 2.725  
f_min = 10e6  
f_max = 10e11 

def planck_spectrum(freq, temp):
    return (8 * pi * h * freq**3) / ((c**3) * (np.exp((h * freq) / (k * temp)) - 1) * (1+z)**3)

redshifts = [0, 6, 20, 1090, 1e4]

#plt.figure(figsize=(12, 8))

for z in redshifts:   
    T = T0 * (1 + z)    
    frequencies = np.linspace(f_min, f_max, 100)    
    intensities = planck_spectrum(frequencies, T)    
    frequencies_ghz = frequencies / 1e9    
    plt.plot(frequencies_ghz, intensities, label=f'Redshift z={z}')

plt.title('Blackbody Radiation Spectrum of the Universe at Different Redshifts')
plt.xlabel('Frequency (GHz)')
plt.ylabel('B')
#plt.yscale('log')
plt.legend()
plt.show()