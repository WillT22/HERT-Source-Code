#%% Import libraries
import numpy as np
import scipy.constants as sc
import scipy.integrate as integrate
import matplotlib.pyplot as plt

#%% Setting up the problem
# Electron mass in MeV:  Conversion from kg*m^2/s^2 (Joules) to MeV.
e_E0 = sc.electron_mass * sc.c**2 / (sc.electron_volt * 1e6)

# Create electron kinetic energy range
KE = np.linspace(0, 2, num = 201)
KE = KE[1:]

# Relativistic momentum of the electron
p = np.sqrt((KE + e_E0)**2 - e_E0**2) / sc.c 

# deBroglie wavelength of an electron
lmbda = (sc.h / (sc.electron_volt * 1e6)) / p 

# Plotting lambda vs KE
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(KE, lmbda)
plt.xlim(0, 2)
#plt.ylim(0, 450)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Electron Wavelength')  # Label the y-axis with the correct unit
plt.title('Electron Wavelength vs Electron Kinetic Energy')
plt.grid(True)
plt.show()

#%% Start with detector placement
#a_D = 18E-3 # (m), detector aperture (HERT)
a_D = 8.382e-3 # (m), detector aperture (Aero)
d_D = 61e-3 + 0.5 # detector opening to last tooth + distance from image slit to detector front

#%% Image Slit
# Theory: what is the smallest opening that allows for the central maximum to fully be detected
theta_I = np.atan(a_D/2 / d_D)
theta_I_deg = theta_I * 180 / np.pi
a_I_min = lmbda / np.sin(theta_I)

# Experimental: what is the maximum order number that can be detected
a_I1 = 15.1e-3
order_I1 = np.pi * a_I1 * np.sin(theta_I) / lmbda / np.pi

# Experimental: what is the maximum order number that can be detected
a_I2 = 5e-3
order_I2 = np.pi * a_I2 * np.sin(theta_I) / lmbda / np.pi

order_I_scale = np.mean(order_I1/order_I2)

# Plotting m vs KE
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(KE, order_I1, label=f"Image Slit = {a_I1 * 10**3:.1f} mm")
plt.plot(KE, order_I2, label=f"Image Slit = {a_I2 * 10**3:.1f} mm")
plt.xlim(0, 2)
plt.yscale('log')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Order Number')
plt.title('Order Number vs Electron Kinetic Energy')
plt.grid(True)
plt.legend()
plt.show()

#%% Object Slit
# How much does the beam spread between the object slit and image slit
d_mag = (507.405 + 299.6 + 605.693) * 10**-3 # Effective length of path through magnetic field
theta_O1 = np.atan(a_I1/2 / d_mag)
theta_O2 = np.atan(a_I2/2 / d_mag)

# Theory: what is the smallest opening that allows for the central maximum to fully be detected
a_O1_min = lmbda / np.sin(theta_O1)
a_O2_min = lmbda / np.sin(theta_O2)

# Experimental: what is the maximum order number that makes it through the image slit
a_O1 = 15.1e-3
a_O2 = 5e-3
order_O1 = np.pi * a_O1 * np.sin(theta_O1) / lmbda / np.pi
order_Omix = np.pi * a_O1 * np.sin(theta_O2) / lmbda / np.pi
order_O2 = np.pi * a_O2 * np.sin(theta_O2) / lmbda / np.pi

order_O_scale12 = np.mean(order_O1/order_O2)
order_O_scale1mix = np.mean(order_O1/order_Omix)

# Plotting m vs KE
plt.figure(figsize=(8, 6))  # Adjust figure size as needed
plt.plot(KE, order_O1, label=f"Object Slit = {a_O1 * 10**3:.1f} mm, Image Slit = {a_I1 * 10**3:.1f} mm")
plt.plot(KE, order_Omix, label=f"Object Slit = {a_O1 * 10**3:.1f} mm, Image Slit = {a_I2 * 10**3:.1f} mm")
plt.plot(KE, order_O2, label=f"Object Slit = {a_O2 * 10**3:.1f} mm, Image Slit = {a_I2 * 10**3:.1f} mm")
plt.xlim(0, 2)
plt.yscale('log')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Order Number')
plt.title('Order Number vs Electron Kinetic Energy')
plt.grid(True)
plt.legend()
plt.show()

#%% Numeric Power integration (unused, computationally too expensive)
#order_range = np.linspace(-np.pi, np.pi, num = 1000)
#sin_func = np.sin(order_range)**2 / order_range**2
#I_O = integrate.trapezoid(np.sin(order_range)**2 / order_range**2, order_range) / np.pi

#%% Energy resolution from magnetic spectrometry
r = 191e-3 # radius of curverture in mm
r_min1 = r - a_I1/2 # minimum radius of curverture that can be seen
r_max1 = r + a_I1/2 # maximum radius of curverture that can be seen

r_min2 = r - a_I2/2 # minimum radius of curverture that can be seen
r_max2 = r + a_I2/2 # maximum radius of curverture that can be seen

# Calculate the magnetic field required to bend an electron beam around a set radius of curvature
B = p * (sc.electron_volt * 1e6) / (sc.elementary_charge * r)

# Reverse that: Calculate the max and min momentum that can be seen by the detector
p_min1 = (B * sc.elementary_charge * r_min1) / (sc.electron_volt * 1e6)
p_max1 = (B * sc.elementary_charge * r_max1) / (sc.electron_volt * 1e6)

p_min2 = (B * sc.elementary_charge * r_min2) / (sc.electron_volt * 1e6)
p_max2 = (B * sc.elementary_charge * r_max2) / (sc.electron_volt * 1e6)

# Find Kinetic Energy from momentum
KE_min1 = np.sqrt((p_min1 * sc.c)**2 + e_E0**2) - e_E0
KE_max1 = np.sqrt((p_max1 * sc.c)**2 + e_E0**2) - e_E0

KE_min2 = np.sqrt((p_min2 * sc.c)**2 + e_E0**2) - e_E0
KE_max2 = np.sqrt((p_max2 * sc.c)**2 + e_E0**2) - e_E0

Energy_resolution1 = (KE_max1 - KE_min1) / KE

Energy_resolution2 = (KE_max2 - KE_min2) / KE

resolution_scale = Energy_resolution1/Energy_resolution2

# Plotting energy resolution vs KE
plt.figure(figsize=(8, 6)) 
plt.plot(KE, Energy_resolution1, label=f"Image Slit = {a_I1 * 10**3:.1f} mm")
plt.plot(KE, Energy_resolution2, label=f"Image Slit = {a_I2 * 10**3:.1f} mm")
plt.xlim(0, 2)
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Energy Resolution')
plt.title('Energy Resolution vs Electron Kinetic Energy')
plt.grid(True)
plt.legend()
plt.show()

#%% Energy resolution as a function of image slit width
r = 191e-3 # radius of curverture in mm

a = np.linspace(0, 20e-3, num = 201) # aperture width from 0 to 20 mm
KE_2d, a_2d = np.meshgrid(KE, a)

r_min = r - a_2d/2 # minimum radius of curverture that can be seen
r_max = r + a_2d/2 # maximum radius of curverture that can be seen

# Calculate the magnetic field required to bend an electron beam around a set radius of curvature
p_2d = np.sqrt((KE_2d + e_E0)**2 - e_E0**2) / sc.c 
B_2d = p_2d * (sc.electron_volt * 1e6) / (sc.elementary_charge * r)

# Reverse that: Calculate the max and min momentum that can be seen by the detector
p_min = (B_2d * sc.elementary_charge * r_min) / (sc.electron_volt * 1e6)
p_max = (B_2d * sc.elementary_charge * r_max) / (sc.electron_volt * 1e6)

# Find Kinetic Energy from momentum
KE_min = np.sqrt((p_min * sc.c)**2 + e_E0**2) - e_E0
KE_max = np.sqrt((p_max * sc.c)**2 + e_E0**2) - e_E0

KE_spread = KE_max - KE_min

Energy_resolution = (KE_max - KE_min) / KE_2d

# Plotting the heatmap
plt.figure(figsize=(8, 6))
plt.imshow(Energy_resolution, extent=[KE.min(), KE.max(), a.min()*10**3, a.max()*10**3],
           aspect='auto', origin='lower', cmap='viridis')  # Use 'viridis' colormap
plt.colorbar(label='Energy Resolution (ΔE/E)')
plt.xlabel('Kinetic Energy (MeV)')
plt.ylabel('Image Slit Width (mm)')
plt.title('Energy Resolution vs. Kinetic Energy and Aperture Width')
plt.show()