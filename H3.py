import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

#E = np.genfromtxt('DMC.csv', delimiter=',')
E6 = np.genfromtxt('DMC_6dim.csv', delimiter=',')
#xs = np.genfromtxt('x.csv', delimiter=',')

"""t = np.linspace(0, len(E[:, 0]), len(E[:, 0]))

plt.plot(t[:2500]*0.02, E[:2500, 2], label='\u0394\u03C4 = 0.02')
plt.plot(t[:2500]*0.02, 200 * np.ones(len(t[:2500])), 'r--', label='N0=200')

plt.xlabel('\u03C4')
plt.ylabel("Number of walkers")
plt.legend(loc='best')
plt.savefig("Walkers_cut.jpg")
plt.show()"""

"""plt.plot(t[:2500]*0.02, E[:2500, 0], label='\u0394\u03C4 = 0.02')
plt.plot(t[:2500]*0.02, 0.375 * np.ones(len(t[:2500])), 'r--', label='Analytical result (E=0.375)')
a = np.std(E[:, 1])
plt.fill_between(t[:2500]*0.02, E[:2500, 0]-a, E[:2500, 0]+a, alpha=0.3, label=r'$\langle E_T \rangle\pm$ one standard deviation')

plt.xlabel('\u03C4')
plt.ylabel(r'$\langle E_T \rangle$ (Ha)')
plt.legend(loc='best')
plt.savefig("Energy_avg_cut.jpg")
plt.show()"""

"""plt.plot(t[:2500]*0.02, E[:2500, 1], label='\u0394\u03C4 = 0.02')
plt.plot(t[:2500]*0.02, 0.375 * np.ones(len(t[:2500])), 'r--', label='Analytical result (E=0.375)')

plt.xlabel('\u03C4')
plt.ylabel(r'$E_T$ (Ha)')
plt.legend(loc='best')
plt.savefig("Energy_cut.jpg")
plt.show()"""

"""nb = 20

a = -5
b = 15
s = (b - a) / nb

N = np.zeros(nb)

for i in range(nb):
    print(i)
    for j in range(7500):
        for x in xs[:, j]:
            if a + i * s < x <= a + (i+1) * s:
                N[i] += 1

phi = N / np.sqrt(np.sum(N**2))

plt.scatter(np.linspace(a, b, nb), phi, label='DMC simulation', marker='o')

x_line = np.linspace(a, b, 200)
phi_an = np.sqrt(2) * np.exp(-np.exp(-x_line) - x_line / 2)

plt.plot(x_line, phi_an, 'k', label='Analytical result')
plt.legend(loc='best')
plt.xlabel('x')
plt.ylabel('\u03A6(x)')
plt.savefig('Wave_function.jpg')

plt.show()"""

plt.figure(figsize=(8, 6))
t = np.linspace(0, len(E6[:, 0]), len(E6[:, 0]))
plt.plot(t[:], E6[:, 1], label='\u0394\u03C4 = 0.01')
plt.plot(t[:], -2.903*np.ones(len(t)), 'r--', label='Analytical result (E=-2.903)')

plt.xlabel('Iterations')
plt.ylabel(r'$\langle E_T \rangle$  (a.u.)')
plt.legend(loc='best')

#plt.savefig('Energy_6dim.jpg')
plt.figure(figsize=(8, 6))
t = np.linspace(0, len(E6[:, 0]), len(E6[:, 0]))
plt.plot(t[:], E6[:, 2], label='\u0394\u03C4 = 0.01')
plt.plot(t[:], 1000*np.ones(len(t)), 'r--', label='N0=1000')

plt.xlabel('Iterations')
plt.ylabel('Number of walkers')
plt.legend(loc='best')

#plt.savefig('Walkers_6dim.jpg')
plt.show()
