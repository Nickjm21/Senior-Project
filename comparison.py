import numpy as np
import time_evolve
import matplotlib.pyplot as plt
timeEvo = time_evolve.timeEvo

hbar = 1
ko = 1
m = 1
tao = (2 * m * (1 ** 2)) / (hbar)

L = 100
N = 1000
dx = L / N
dt = 0.1

wavfuncinitial = np.zeros((N, 1), dtype=complex)
wavfuncfinal = np.zeros((N, 1), dtype=complex)
diffarray = np.zeros((N, 1), dtype=complex)


def initial(a, x):
    val = (a ** (-1/2)) * ((2*np.pi) ** (-1/4)) \
        * (np.exp(1j * ko * x)) \
        * (np.exp((-x ** 2) / (4 * (a ** 2))))
    return val


def actualEvo(a, t, T, x):
    v = (a ** (-1/2)) * ((2 * np.pi) ** (-1/4)) * ((1 + (1j * t)/T) ** (-1/2))\
        * (np.exp(1j * (T/t) * ((x/2*a) ** 2))) \
        * (np.exp((((-1j * T)/(4 * (a ** 2) * t))
        * ((x - hbar * ko * t) / m) ** 2) \
        / ((1 + 1j * t) / T)))
    return v


for i in range(N):
    wavfuncinitial[i, :] = initial(1, (i-(N/2)) * dx)


def finalarray(t):
    for i in range(N):
        wavfuncfinal[i, :] = actualEvo(1, t, tao, (i-(N/2)) * dx)
    return wavfuncfinal


def difference(t):
    for i in range(N):
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
