import matplotlib.pyplot as plt
import numpy as np

E = np.genfromtxt('DMC.csv')

t = np.linspace(0, len(E), len(E))

plt.plot(t[11000:], E[11000:])
print(E[11000:].mean())
plt.show()