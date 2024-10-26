import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Sample x-y dataset
#x = np.array([1, 2, 3, 4, 5])
#y = np.array([2, 4, 1, 5, 3])

x, y = [], []
with open('XX.dat', 'r') as f1:
    for line in f1:
        data = line.strip().split()
        x.append(float(data[0]))
        y.append(float(data[1]))

# Create spline interpolation function
spline_func = interp1d(x, y, kind='cubic')

# Generate additional points between each pair of data points
new_x = np.linspace(x[0], x[-1], num=20001)  # Adjust the 'num' parameter for desired density of points

# Use the spline function to obtain corresponding y-values
new_y = spline_func(new_x)

data = np.column_stack((new_x, new_y))
np.savetxt('XX_spline.dat', data, fmt='%.6f', delimiter='\t')

# Print the generated data
#for i in range(len(new_x)):
#    print(f'({new_x[i]:.2f}, {new_y[i]:.2f})')

plt.plot(x,y)
plt.plot(new_x,new_y)
plt.xlabel("time (fs)")
plt.ylabel("nu_d")
plt.savefig("XX_data.png")
