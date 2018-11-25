import numpy as np
import time_evolve
import matplotlib.pyplot as plt
import config as c
timeEvo = time_evolve.timeEvo

# Variables:
tao = (2 * c.mass * (1 ** 2)) / (c.hbar)


wavfuncinitial = np.zeros((c.N, 1), dtype=complex)
wavfuncfinal = np.zeros((c.N, 1), dtype=complex)
diffarray = np.zeros((c.N, 1), dtype=complex)


def initial(a, x):
    val = (a ** (-1/2)) * ((2*np.pi) ** (-1/4)) \
        * (np.exp(1j * c.ko * x)) \
        * (np.exp((-x ** 2) / (4 * (a ** 2))))
    return val


def actualEvo(a, t, T, x):
    v = (a ** (-1/2)) * ((2 * np.pi) ** (-1/4)) * ((1 + (1j * t)/T) ** (-1/2))\
        * (np.exp(1j * (T/t) * ((x/2*a) ** 2))) \
        * (np.exp((((-1j * T)/(4 * (a ** 2) * t))
            *((x - c.hbar * c.ko * t) / c.mass) ** 2)
            / ((1 + 1j * t) / T)))
    return v


for i in range(c.N):
    wavfuncinitial[i, :] = initial(1, (i-(c.N/2)) * c.dx)


def finalarray(t):
    for i in range(c.N):
        wavfuncfinal[i, :] = actualEvo(1, t, tao, (i-(c.N/2)) * c.dx)
    return wavfuncfinal


def difference(t):
    for i in range(c.N):
        diffarray[i, :] = finalarray(t)[i, :] \
            - timeEvo(t, wavfuncinitial)[i, :]
    return diffarray


plt.plot(abs(difference(1)))

plt.plot(abs(timeEvo(100, wavfuncinitial)))
plt.plot(abs(finalarray(10)))

plt.plot(abs(timeEvo(10, wavfuncinitial)))
plt.plot(abs(timeEvo(100, wavfuncinitial)))
plt.plot(abs(timeEvo(150, wavfuncinitial)))
plt.plot(abs(timeEvo(200, wavfuncinitial)))
plt.plot(abs(finalarray(1)))
plt.plot(abs(finalarray(10)))
plt.plot(abs(finalarray(15)))
plt.plot(abs(finalarray(20)))
