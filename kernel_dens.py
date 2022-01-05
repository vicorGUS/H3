import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import KernelDensity

"""E_all = np.genfromtxt('DMC.csv', delimiter=',')
i_start = 11000
ET = E_all[i_start:, 0]
ET_avg = E_all[i_start:, 1]
ET_avg1 = E_all[i_start:, 2]
t = np.linspace(0, len(ET), len(ET))
print(ET_avg[-1])"""

x_raw = np.genfromtxt('x.csv')
x = np.array(sorted(x_raw, key = lambda a:float(a)))[:, np.newaxis]
x_plot = np.linspace(min(x), max(x), len(x))

phi = np.sqrt(2) * np.exp(-np.exp(-x) - x/2)
kde = KernelDensity(kernel="gaussian", bandwidth=(i+1)/10).fit(x)
log_dens = kde.score_samples(x_plot)

plt.figure()
plt.plot(x, phi**2., color='k', label='Phi^2')
plt.plot(x, phi, color='r', label='Phi')
plt.hist(x, bins=20, density=True, color='r', alpha=0.5, label='Histogram')
plt.plot(x_plot, np.exp(log_dens), color='b', label='Kernel density')
plt.legend()
