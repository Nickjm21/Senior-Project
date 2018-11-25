import numpy as np
import cmd
import matplotlib.pyplot as plt
import potentials
import abc
import config as c
from TDMsolver import TDMAsolver
TDMA = TDMAsolver

'''
This file is a Crank-Nicolson (CN) approximation to the time evolution
operator in QM. Units are as follows:
- Length in nanometers (nm)
- Time in femtoseconds (fm)
- Energies in electronvolts (eV)
'''

abc.num

print(potentials.potfuncs)
print(TDMAsolver)
# Establishing potential function
pot = 0

# Elements of the time evo matrix:

alpha = (c.dt * c.hbar)/(4 * c.mass * (c.dx**2))

ui = (c.dt / (2 * c.hbar)) * pot


# Functions that the initial wavefunction can take
def Gaussian(x):
    val = np.exp(-(x ** 2))
    return(val)


def NormGauss(a, x):
    val = (a ** (-1/2)) * ((2*np.pi) ** (-1/4)) \
        * (np.exp(1j * c.ko * x)) \
        * (np.exp((-x ** 2) / (4 * (a ** 2))))
    return val


def MMGauss(x):
    val = ((2 * np.pi * (c.deltax ** 2)) ** (-1/4)) \
        * (np.exp(
            -((x - c.xo) ** 2) / (4 * (c.deltax ** 2))
            + 1j * c.ko * x
        ))
    return val


# Checking normalization function
def checknorm(func):
    norm = 0.
    for i in range(c.N):
        norm = (abs(func[i, :]) ** 2) + norm
    value = 1 - (c.dx * norm)
    return value


# Allocating full matricies
tevo1 = np.zeros((c.N, c.N), dtype=complex)
tevo2 = np.zeros((c.N, c. N), dtype=complex)
wavfunc = np.zeros((c.N, 1), dtype=complex)

# Allocating parts of matricies for TDMA
tevo1top = np.zeros((c.N-1, 1), dtype=complex)
tevo1diag = np.zeros((c.N, 1), dtype=complex)
tevo1bot = np.zeros((c.N-1, 1), dtype=complex)
# tevo2top = np.zeros((N-1, 1), dtype=complex)
# tevo2diag = np.zeros((N, 1), dtype=complex)
# tevo2bot = np.zeros((N-1, 1), dtype=complex)

# Filling tevo2
for i in range(c.N):
    if i == 0:
        tevo2[i, :] = [
            1 - 2j * alpha - 1j * ui if j == 0
            else 1j * alpha if j == 1
            else 0 for j in range(c.N)]
    if i == c.N:
        tevo2[i, :] = [
            1j * alpha if j == c.N - 1
            else 1 - 2j * alpha - 1j * ui if j == c.N
            else 0 for j in range(c.N)]
    else:
        tevo2[i, :] = [
            1j * alpha if j == i - 1 or j == i + 1
            else 1 - 2j * alpha - 1j * ui if j == i
            else 0 for j in range(c.N)]

print(tevo2)

# Filling the wavefunction
for i in range(c.N):
    wavfunc[i, :] = MMGauss((i-(c.N / 2)) * c.dx)

# Filling in the top and bottom parts of the two matricies for the TDMA solver.
# tevo2 doesn't need to be filled for this part because tevo1 is the matrix
# that will be used in the TDMA solver.
for i in range(c.N-1):

    tevo1bot[i, :] = -alpha*1j
    tevo1top[i, :] = -alpha*1j
    # tevo2bot[i, :] = alpha*1j
    # tevo2top[i, :] = alpha*1j

# Filling in the middle diagonal parts of the two matricies for the TDMA solver
# Again tevo2 doesn't need to be filled this way.
for i in range(c.N):
    tevo1diag[i, :] = 1 + 2j*alpha + 1j*ui
    # tevo2diag[i, :] = 1 - 2j*alpha - 1j*ui


def timeEvo(num, wavefunction):
    for i in range(num):
        partialevo = np.matmul(tevo2, wavefunction)
        evowavefunc = TDMA(tevo1bot, tevo1diag, tevo1top, partialevo)
        wavefunction = evowavefunc

    return evowavefunc


timeEvo(10, wavfunc)

plt.plot(abs(wavfunc))
plt.plot(abs(timeEvo(100, wavfunc)))

for i in range(1, 10):
    if np.greater(checknorm(timeEvo(i, wavfunc))[0], np.float64(10. ** -13)):
        print(i)
