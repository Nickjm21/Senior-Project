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


for i in range(c.N):
    wavfunc[i, :] = TE.MMGauss((i-(c.N / 2)) * c.dx, c.ko)


time300 = TE.timeEvo(300, wavfunc)

time3002 = TE.timeEvo(300, wavfunc)


def probtrans(wavefunction):
    norm = 0.
    for i in range(c.N):
        if i > ((c.potmax / c.dx) + c.N / 2):
            norm = (abs(wavefunction[i, :]) ** 2) + norm
    value = (c.dx * norm)[0]
    return value


probtrans(time300)
probtrans(time3002)


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


Earray = np.arange(0.0, 5, 0.001) / c.Vo
T = transmition(Earray)

fig, ax = plt.subplots()
ax.plot(Earray, T)

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')
ax.grid()

fig.savefig("test.png")
plt.show()


probtrans(abs((TE.timeEvo(300, TE.wavfunc))))


initialfuncarray = {}
for i in range(4):
    initialfuncarray[i] = FillWavfunc(i)


# Expiremental transmission probability:
krange = np.arange(0.1, 8, 0.1)
Erange = Energy(krange)

krange
Erange

initialfuncarray = {}
for i in krange:
    wavfunc = np.zeros((c.N, 1), dtype=complex)
    for j in range(c.N):
        wavfunc[j, :] = TE.MMGauss((j-(c.N / 2)) * c.dx, i * c.ko)
        initialfuncarray[int(round(i / 0.1))] = wavfunc

timedarray = {}
for i in krange:
    print(int(round(i / 0.1)))
    timedarray[int(round(i / 0.1))] = TE.timeEvo(300, initialfuncarray[int(round(i / 0.1))])

transarray = {}
for i in krange:
    transarray[int(round(i / 0.1))] = probtrans(timedarray[int(round(i / 0.1))])

transarray[1]

Tarray = np.zeros((79, 1), dtype=float)
for i in range(1, 79):
    Tarray[i, :] = transarray[i]


# @@@@@@@@@@ PLOTS @@@@@@@@@@


ax.plot()
plt.plot(Tarray)

fig, ax = plt.subplots()

plt.plot(Erange, Tarray)
# plt.errorbar(krange, Tarray, c.deltaE(krange))

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')
ax.grid()
plt.xlim(0, 0.5)

fig.savefig("test.png")
plt.show()


ax.plot()
plt.plot(Tarray)

fig, ax = plt.subplots()

plt.plot(10 * Erange, Tarray)
# plt.errorbar(krange, Tarray, c.deltaE(krange))

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')
ax.grid()

fig.savefig("test.png")
plt.show()



Earray = np.arange(0.0, 5, 0.001)
T = transmition(Earray)

plt.errorbar(Earray, T * 10, 0.1)
ax.plot(Earray, T)

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')

ax.grid()

fig.savefig("test.png")
plt.show()


fig = plt.figure(0)
x = np.arange(10.0)
y = np.sin(np.arange(10.0) / 20.0 * np.pi)

upperlimits = np.array([1, 0] * 5)
lowerlimits = np.array([0, 1] * 5)

plt.xlim(-1, 10)


x = krange
y = Tarray

plt.scatter(x, y, 1.7, "orange")
plt.xlabel("Leprechauns")
plt.ylabel("Gold")
plt.show()





Earray = np.arange(0.0, 10, 0.001) / c.Vo
T = transmition(Earray)

fig, ax = plt.subplots()
ax.plot(Earray, T)

ax.set(xlabel='E / Vo', ylabel='Transmition Probability',
       title='Transmition through a square barrier')
ax.grid()

fig.savefig("test.png")
plt.show()



plt.scatter(krange * 1.7, y, 1.7, "orange")
plt.plot(Earray, T)

plt.xlim(0, 4)

plt.show()
