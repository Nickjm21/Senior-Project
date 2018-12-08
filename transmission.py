import numpy as np
import matplotlib.pyplot as plt
import cmath
import config as c
import time_evolve as TE


def newsqrt(x):
    return cmath.sqrt(x)


sqrt2 = np.vectorize(newsqrt)

# Initializing wavefunctions:

wavfunc = np.zeros((c.N, 1), dtype=complex)


def FillWavfunc(k):
    for i in range(c.N):
        wavfunc[i, :] = TE.MMGauss((i-(c.N / 2)) * c.dx, k)
    return wavfunc


def probtrans(wavefunction):
    norm = 0.
    for i in range(c.N):
        if i > ((c.potmax / c.dx) + c.N / 2):
            norm = (abs(wavefunction[i, :]) ** 2) + norm
    value = (c.dx * norm)[0]
    return value


def probreflect(wavefunction):
    norm = 0.
    for i in range(c.N):
        if i < ((c.potmin / c.dx) + c.N / 2):
            norm = (abs(wavefunction[i, :]) ** 2) + norm
    value = (c.dx * norm)[0]
    return value


def probinpot(wavefunction):
    norm = 0.
    for i in range(c.N):
        if i > ((c.potmax / c.dx) + c.N / 2):
            norm = (abs(wavefunction[i, :]) ** 2) + norm
    value = (c.dx * norm)[0]
    return value


def Energy(k):
    val = ((c.hbar ** 2) * ((k ** 2) - (1 / (4 * (c.deltax ** 2))))) \
        / (2 * c.mass)
    return val


# Theory transmition probability

def transmition(E):
    val = 1 / (1 + ((c.Vo ** 2)
            * (np.sinh(c.a * sqrt2(2 * c.mass * (c.Vo - E) / c.hbar))) ** 2) / (
        4 * E * (c.Vo - E)))
    return np.real(val)


Earray = np.arange(0.0, 4, 0.001)
T = transmition(Earray)

fig, ax = plt.subplots()
ax.plot(Earray, T)

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')
ax.grid()

fig.savefig("test.png")
plt.show()


probtrans(abs((TE.timeEvo(300, TE.wavfunc))))


# Expiremental transmission probability:

initialfuncarray = {}
for i in range(4):
    initialfuncarray[i] = FillWavfunc(i)

timedarray = {}
for i in range(4):
    timedarray[i] = TE.timeEvo(300, initialfuncarray[i])

transarray = {}
for i in range(4):
    transarray[i] = probtrans(timedarray[i])

transarray
