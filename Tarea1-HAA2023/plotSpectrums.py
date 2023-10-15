# Script Plot Sa(T)

# Import
import matplotlib.pyplot as plt

# Data
T = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
Sx = [0.110, 0.199, 0.294, 0.304, 0.260, 0.210, 0.168, 0.137, 0.114, 0.096, 0.082, 0.072, 0.063, 0.056, 0.050, 0.045, 0.041, 0.037, 0.034, 0.031, 0.029]
Sy = [0.096, 0.172, 0.254, 0.263, 0.225, 0.181, 0.146, 0.119, 0.099, 0.083, 0.071, 0.062, 0.054, 0.048, 0.043, 0.039, 0.035, 0.032, 0.030, 0.027, 0.025]

# Plot
plt.figure(figsize=(10, 5))
plt.plot(T, Sx, marker='o', linestyle='-', label='Sx(T)')
plt.plot(T, Sy, marker='o', linestyle='-', label='Sy(T)')
plt.title('Sa(T)')
plt.xlabel('T [s]')
plt.ylabel('Sa(T)_reducido/g')
plt.grid(True)
plt.legend()
plt.show()
