import numpy as np
from scipy.optimize import least_squares

## Importing data from MATLAB ##
count_rate = np.loadtxt("E:/HERT_Drive/Matlab Main/Bow Tie/count_rate.txt")
geo_factor = np.loadtxt("E:/HERT_Drive/Matlab Main/Bow Tie/geo_factor.txt")
bin_midpoints = np.loadtxt("E:/HERT_Drive/Matlab Main/Bow Tie/energy_bins.txt")
bowtie_flux = np.loadtxt("E:/HERT_Drive/Matlab Main/Bow Tie/bowtie_flux.txt")

## Defining least-squares function integral ##
def func(flux,count_rate,geo_factor):
    return np.trapz(flux * geo_factor, bin_midpoints)
            
def chisq(flux, count_rate, geo_factor):
    integral = func(flux,count_rate,geo_factor)
    return np.subtract(count_rate, integral)

channel_flux = np.empty(len(count_rate))
for channel in range(len(count_rate)):
   print('Channel Number', channel+1)
   channel_flux_temp = least_squares(chisq, bowtie_flux[channel], args=(count_rate[channel], geo_factor[channel]), method='lm')
   channel_flux[channel] = channel_flux_temp.x
   
np.savetxt("Bow Tie/lsqr_flux.txt",channel_flux)