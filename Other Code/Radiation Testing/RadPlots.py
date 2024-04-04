#!/usr/bin/env python
# coding: utf-8

# # Data Analysis of 0-65 krad ADC Test

# In[1]:


# import modules
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import csv
import matplotlib.cm as cm
cmap = cm.get_cmap('plasma')


# In[2]:


# Specify the path to your text file
file_path = "RadTestResults_total.txt"

# Read the data into a DataFrame
data = pd.read_csv(file_path, sep=',',header=None)
# Set the row and column names to their default integer index
data.index = range(len(data))
data.columns = range(len(data.columns))

# Filter the DataFrame to remove rows containing the specified string in any column
data = data[~data[0].str.contains("Error detected at")]
data[0] = pd.to_numeric(data[0], errors='coerce')


# In[3]:


numpy_array = data.to_numpy()

# Extract the first 8 columns as 'stddev'
stddev = numpy_array[:, :8]

# Extract the next 8 columns as 'means'
means = numpy_array[:, 8:16]

# Extract the last column as time data in a vector
time_data = numpy_array[:, -1]

# Convert time strings to datetime objects
time_data = [datetime.strptime(time_str.strip(), "%d-%b-%Y %H:%M:%S") for time_str in time_data]


# ### Calculate CI for Radiation Dosage of Duration

# In[5]:


# Copy CI of radiation rate based on distance from Radiation_Calibration.ipynb
def ci_rate(distance: np.ndarray):
    # Calculate radiation rate using 'a' and 'b'
    rate = a * (distance ** b)

    # Check if the input 'distance' is a single float or int
    if isinstance(distance, (float, int)):
        confidence_intervals = np.zeros((1, 3))
        if distance < 1:
            # Calculate lower and upper bounds of radiation rate when distance is less than 1
            confidence_intervals[0, 0] = conf_int_a[0] * (distance ** conf_int_b[1])
            confidence_intervals[0, 2] = conf_int_a[1] * (distance ** conf_int_b[0])
        else:
            # Calculate lower and upper bounds of radiation rate when distance is greater than or equal to 1
            confidence_intervals[0, 0] = conf_int_a[0] * (distance ** conf_int_b[0])
            confidence_intervals[0, 2] = conf_int_a[1] * (distance ** conf_int_b[1])
    
    # When 'distance' is an array
    else:
        # initialize a 2D NumPy array for confidence intervals
        confidence_intervals = np.zeros((len(distance), 3))
        # Propagate uncertainty through the power equation for each distance value in the input array
        for i in range(len(distance)):
            # The middle column (1) contains the predicted radiation rate
            confidence_intervals[i, 1] = rate[i]
            if distance[i] < 1:
                # Calculate lower and upper bounds when distance is less than 1
                confidence_intervals[i, 0] = conf_int_a[0] * (distance[i] ** conf_int_b[1])
                confidence_intervals[i, 2] = conf_int_a[1] * (distance[i] ** conf_int_b[0])
            else:
                # Calculate lower and upper bounds when distance is greater than or equal to 1
                confidence_intervals[i, 0] = conf_int_a[0] * (distance[i] ** conf_int_b[0])
                confidence_intervals[i, 2] = conf_int_a[1] * (distance[i] ** conf_int_b[1])

    return confidence_intervals


# In[12]:


# Define the reference datetime (11:59 AM on October 30, 2023)
reference_datetime = datetime(2023, 10, 30, 11, 59)

# Calculate time differences in seconds
time_in_seconds = np.array([(dt - reference_datetime).total_seconds() for dt in time_data])

# Importing Coefficients and Confidence interval from R
a = 0.1433639
b = -1.991746 

conf_int_a = np.array([0.1383645, 0.1483145])
conf_int_b = np.array([-2.044445, -1.941195])

# Determining max and min rates
rate = np.array([ci_rate(0.42)[0,0], a * (0.4 ** b), ci_rate(0.38)[0,2]])

# Determine the expected radiation dose at each time
rad_dose = rate[1] * time_in_seconds
# Determining max and min possible radiation dosages
rad_dose_end = rate * time_in_seconds[-1]

print('Minimum Radiation Dose = %.2f rad' % rad_dose_end[0])
print('Expected Radiation Dose = %.2f rad' % rad_dose_end[1])
print('Maximum Radiation Dose = %.2f rad' % rad_dose_end[2])


# In[13]:


# Create a figure for the standard deviation data
fig1, ax1 = plt.subplots(figsize=(16, 8))
ax1.set_ylabel('Standard Deviation (DN)')
for channel in range(8):
    color = cmap(channel / 8)[0:7]  # Assign a color based on the channel (0 to 7)
    ax1.plot(time_data, stddev[:, channel], label=f'Channel {channel + 1}', color=color)
ax1.legend(loc='upper left')
ax1.set_xlabel('Time')

# Create a second X-axis for radiation data
ax2 = ax1.twiny()
ax2.set_xlabel('Expected Radiation Exposure (rad)')
# Plot the radiation data on the top X-axis
ax2.plot(rad_dose, [stddev[0, 0]] * len(rad_dose), alpha=0)  # Create a hidden line for the top X-axis

fig1.autofmt_xdate()

# Set the title and add some space between the title and the radiation axis label
fig1.suptitle('Standard Deviation Over Time', y=1.0)

plt.show()


# In[11]:


# Create a figure for the standard deviation data
fig2, ax3 = plt.subplots(figsize=(16, 8))
ax3.set_ylabel('Mean (DN)')
for channel in range(8):
    color = cmap(channel / 8)[0:7]  # Assign a color based on the channel (0 to 7)
    ax3.plot(time_data, means[:, channel], label=f'Channel {channel + 1}', color=color)
ax3.legend(loc='upper left')
ax3.set_xlabel('Time')

# Create a second X-axis for radiation data
ax4 = ax3.twiny()
ax4.set_xlabel('Expected Radiation Exposure (rad)')
# Plot the radiation data on the top X-axis
ax4.plot(rad_dose, [means[0, 0]] * len(rad_dose), alpha=0)  # Create a hidden line for the top X-axis

fig2.autofmt_xdate()

# Set the title and add some space between the title and the radiation axis label
fig2.suptitle('Mean Over Time', y=1.0)

plt.show()
