import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# fwhm=0.15
# fwhm / (2 * np.sqrt(2 * np.log(2)))


def sigma(x, par1, par2):
    ''' From J. Osán, J. de Hoog, P. Van Espen, I. Szalóki, C. ‐U. Ro, R. Van Grieken, Evaluation of energy‐dispersive x‐ray spectra of low‐ Z elements from electron‐probe microanalysis of individual particles, X-Ray Spectrom. 30 (2001) 419–426. https://doi.org/10.1002/xrs.523.'''
    sigma_val = np.sqrt(par1 + par2 * x)
    return sigma_val

def sigma_from_fwhm(fwhm):
    sigma = fwhm/(2 * np.sqrt(2 * np.log(2)))
    return sigma

# Given values
C_sigma = 2.20e-02
C_x = 0.2774


# All peaks with energy lower than 3keV have a tail
data = {'B_Ka1': {'sigma' : 0.0185, 'energy' : 0.182},
        'O_Ka1': {'sigma' : 0.0235, 'energy' : 0.523},
        'Mn_Ka1': {'sigma' : sigma_from_fwhm(0.126), 'energy' : 5.8987},
        'Ti_Ka1': {'sigma' : 0.048, 'energy' : 4.507},
        'P_Ka1': {'sigma' : 0.035, 'energy' : 2.006},
        'Co_Ka1': {'sigma' : 0.058, 'energy' : 6.921},
        'Na_Ka1': {'sigma' : 0.028, 'energy' : 1.034},
        'Cl_Ka1': {'sigma' : 0.038, 'energy' : 2.620},
        'In_La1': {'sigma' : 0.043, 'energy' : 3.28}
        }



data_sigma = [value['sigma'] for value in data.values()]
data_x = [value['energy'] for value in data.values()]



# Fit the sigma function to the data points
initial_guess = [1e-6, 1e-6]  # Initial guess for par1 and par2
popt, pcov = curve_fit(sigma, data_x, data_sigma, p0=initial_guess)

# Extract fitted parameters
par1, par2 = popt
print(f"Fitted par1: {par1}")
print(f"Fitted par2: {par2}")



# Generate x values for the plot (for smoother curve)
x_fit = np.linspace(min(data_x), max(data_x), 100)
sigma_fit = sigma(x_fit, *popt)

# Plot the original data points and the fitted curve
plt.figure(figsize=(8, 6))
plt.plot(data_x, data_sigma, 'bo', label='Data points')  # Original data points
plt.plot(x_fit, sigma_fit, 'r-', label='Fitted curve')   # Fitted curve
for key, value in data.items():
    plt.text(value['energy'], value['sigma']+0.001, key)
plt.xlabel('Energy (x)')
plt.ylabel('Sigma')
plt.title('Fitting Sigma vs Energy')
plt.legend()
plt.grid(True)
plt.show()

# Find detector params
fano_factor = 0.12 # From Ritchie, 2009
conv_eff = par2 /fano_factor
elec_noise = np.sqrt(par1) / conv_eff
print(f"Electronic noise : {elec_noise:.5f}")
print(f"Conversion efficiency : {conv_eff*1000:.2f} eV") #Matches perfectly the result from Ritchie!
print(f"Conversion efficiency : {conv_eff:.8f} keV") #Matches perfectly the result from Ritchie!