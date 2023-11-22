'''
Created by Yi Jie (Josh) Tseng on 9/28/2023

The purpose of this code is to graph scatter plot and the slop of selected region. This script will take a user input file and 
plot all date into a scatter plot using a column named "Elapsed" as the X-axis and the sample cell counts as the Y-axis. After 
checking the initial scatter plot, you will be asked to input a lower time value and a higher time value. This code will then 
plot the slope lines of the data within the inputted time values. The figures will be saved as {filename.split('.')[0]}_Combined_Scatter.png
and {filename.split('.')[0]}_Regression_lines.png; they should be where the data file is.

To run this script, navigate to the folder where the script is and do:
python3 slope_graphing.py <path/to/input/file>

'''

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress
from sys import argv

# User input the file
filename = ' '.join(argv[1:])
fileN_no_ext = filename.split('.')
fileN_no_ext = '.'.join(fileN_no_ext[:-1])

# Loading data
df = pd.read_csv(filename, delimiter='\t', skiprows=6)

# Get sample names
sample_names = list(df.columns)[2:]

# Create a figure for the combined plot
plt.figure(figsize=(12, 6))

# Create a loop to plot each sample
for sample in sample_names:
    plt.scatter(df['Elapsed'], df[sample], label=sample)

# Add labels and a legend for the combined plot
plt.xlabel('Time (hr)',  fontsize=14)
plt.ylabel('Cell count',  fontsize=14)
plt.title('Combined Scatter Plot',  fontsize=25)
plt.legend(loc='center left', 
           bbox_to_anchor=(1.02, 0.5))
plt.tight_layout()
print(filename.split(".")[0])
plt.savefig(f'{fileN_no_ext}_Combined_Scatter.png', dpi=300)

print('[+] Scatter plot for all data points is generated')
print('[+] Plotting local region...')
print(f'[+] {fileN_no_ext}_Combined_Scatter.png is genereated')


time_min = float(np.min(df.Elapsed))
time_max = float(np.max(df.Elapsed))

# get user inputs
try:
    # Getting user input for lower time value
    min_time = float(input(f'[+] Please input a lower time value ranging from {time_min} to {time_max}:').strip())
    if min_time in df['Elapsed'].values and min_time >= time_min:
        low_index = df[df['Elapsed'] == min_time].index.tolist()[0]
        closest_value_low = df['Elapsed'].iloc[low_index]
    else:
        low_index = df['Elapsed'].sub(min_time).abs().idxmin()
        closest_value_low = df['Elapsed'].iloc[low_index]
        print(f'[+] User input: {min_time} is not found in the data, rounding {min_time} to the nearest integer: {closest_value_low}')
    
    # Getting user input for high time value
    max_time = float(input(f'[+] Please input a higher time value ranging from {min_time} to {time_max}:').strip())
    if max_time in df['Elapsed'].values and max_time <= time_max:
        high_index = df[df['Elapsed'] == max_time].index.tolist()[0]
        closest_value_high = df['Elapsed'].iloc[high_index]
    else:
        high_index = df['Elapsed'].sub(max_time).abs().idxmin()
        closest_value_high = df['Elapsed'].iloc[high_index]
        print(f'[+] User input: {max_time} is not found in the data, rounding {max_time} to the nearest integer: {closest_value_high}')
except ValueError:
    print('[+] Please input a valid value')


df_new = df.iloc[low_index:high_index+1]

# Create a figure for the combined plot
#plt.figure(figsize=(12, 6))

# Create a loop to plot each sample
for sample in sample_names:

    # Calculate linear regression statistics
    result = linregress(df_new['Elapsed'], df_new[sample])

    # Access the slope from the result
    slope = result.slope

    # Create a line using the calculated slope and the x-values
    x_values = np.array(df_new.Elapsed)
    y_values = slope * x_values + result.intercept

    # Plot the regression line
    plt.plot(x_values, y_values, label=sample)


# Set the x-axis limits
x_lower_limit = closest_value_low - 1  # Adjust as needed
x_upper_limit = closest_value_high + 1  # Adjust as needed
plt.xlim(x_lower_limit, x_upper_limit)

# Add labels and legend
plt.xlabel('Time (hr)', fontsize=14)
plt.ylabel('Cell count', fontsize=14)
plt.title(f'Regression Line from {closest_value_low} to {closest_value_high}', fontsize=25)
#plt.legend()
#plt.legend(loc='center left', 
#           bbox_to_anchor=(1.02, 0.5))
plt.tight_layout()

plt.savefig(f'{fileN_no_ext}_Regression_lines.png', dpi=300)

print(f'[+] {fileN_no_ext}_Regression_lines.png is generated')
print('-----------\n[+] Done')
