#!/usr/bin/env python
# coding: utf-8

# # Radiation Callibration Curve and Distances (with error)

# Using equation $y = a x^{b}$, we try to best-fit our data to minimize error and get an accurate reading of radiation as a measure of distance.

# In[1]:


# import modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import statistics as s
import csv


# In[2]:


# input data
distance  = [0.5,1,2,3]
radiation = [0.57031,0.14243,0.037307,0.01711]


# In[3]:


# plot data
dist0 = np.asarray(distance)
rad0  = np.asarray(radiation)
plt.plot(dist0, rad0, 'o')
plt.xlabel('distance (m)')
plt.ylabel('radiation rate (rads/s)')
plt.title("Measured Radiation Data")


# Taken from https://stackoverflow.com/questions/72204658/power-law-fit-in-python to utilize least squares function to find best fit.

# In[4]:


# utilize the least-squares method to find a power equation of best fit in the form y=a*x^b
def fitpower( y, x: np.ndarray, c0=1, p0=1, verbose=0, **kw
) -> "OptimizeResult with .x .fun ...":
    """ fit y ~ c x^p with scipy least_squares """
    # minimize the error
    def err( cp ):
        c, p = cp
        return y - c * (x ** p)
    # calculate residuals
    res = least_squares( err, x0=[c0, p0], jac='3-point', verbose=verbose, **kw )
    c, p = res.x
    if verbose:
        print( "fitpower: y ~ %.3g x^%.3g \n  err %s " % (
            c, p, res.fun ))
    return res  # an OptimizeResult with .x .fun ...


# In[5]:


# use best fit equation to find coeffecients a and b
fit = fitpower(radiation, distance, verbose=1)
a = fit.x[0]
b = fit.x[1]
print("a = %a" % a)
print("b = %a" % b)


# Now, we create a function that takes radiation values and outputs distances from the source:
# 
# $x = (y/a)^{1/b}$

# In[6]:


# create a function to translate radiation rate into distance
def rad(rate: np.ndarray, a, b):
    dist = np.power(rate/a,1/b)
    return dist


# In[7]:


# show calculated distances versus measured distances
rad_test = rad(radiation,a,b)
print(rad_test)
print(distance)


# In[8]:


# find difference between actual and calculated distances
dif = distance-rad_test
# find mean and standard deviation of distances
meandif = s.mean(dif)
stddif  = s.stdev(dif)
# calculate error for use in error bars
rad_err = dif + 2*stddif
print(dif)
print(meandif)
print(stddif)


# In[19]:


# create a smooth distribution of radiation rate data points
rad_in = np.linspace(5,0.01,num=1000)
# and calculate their corresponding distances
dist_in = rad(rad_in,a,b)

# the closest we can get is 0.35m
closest = 0.35
closest_rate = a*np.power(closest,b)
print('Closest rate = ', closest_rate,'rad/s')
rad_in_new = [x for x in rad_in if x < closest]
dist_in_new = rad(rad_in_new,a,b)

# plot data
fig, ax = plt.subplots()
ax.plot(dist_in_new, rad_in_new,label='Predicted Radiation Rate v Distance')
# and plot error bars
ax.errorbar(np.power(rad0/a,1/b), rad0, xerr=rad_err,marker = 'o', ls='', label= 'Prediction of Measured Values')
ax.set(xlabel='distance (m)', ylabel='dose rate (rad/s)')
ax.legend()
plt.title("Radiation Rate v. Distance")
plt.show()


# In[21]:


# test out a value for radiation rate
rate = np.linspace(closest_rate,0.5,50)
dist1 = rad(rate,a,b) 
#print("Distance for %a rads/s: %a" % (rate, rad1))

# find total time to get 100krads with this rate
time = 100e3/rate
#print("time in seconds: %.2f" % time)
time_min = time/60
#print("time in minutes: %.2f" % time_min)

array = np.transpose([[*rate],[*dist1],[*time]])
print(array)


# ### write to csv file
# f = open('radiation_rate.csv', 'w',newline='')
# writer = csv.writer(f)
# writer.writerow(['Rad Rate','Distance','Time'])
# writer.writerows(array)
# f.close()

# In[23]:


with open('radiation_rate_new.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Rad Rate', 'Distance', 'Time'])

    # Write the array to the CSV file with scientific notation
    for row in array:
        formatted_row = [f'{value:.8e}' for value in row]
        writer.writerow(formatted_row)


# In[ ]:




