import numpy as np
import matplotlib.pyplot as plt
import TDMsolver
TDMA = TDMsolver.TDMAsolver

# Defining variables
mass = 1
hbar = 1
ko = 1

L = 100  # Range of x values
N = 1000  # Number of points
dx = L / N  # Size of position steps
nsteps = 10  # Number of time steps
dt = 0.1  # Size of time step
pot = 0  # Potential function (Unsure how to implement an actual function)
alpha = (dt * hbar)/(4 * mass * (dx**2))
ui = (dt / (2 * hbar)) * pot


# Functions that the initial wavefunction can take
def Gaussian(x):
    val = np.exp(-(x ** 2))
    return(val)


def NormGauss(a, x):
    val = (a ** (-1/2)) * ((2*np.pi) ** (-1/4)) \
        * (np.exp(1j * ko * x)) \
        * (np.exp((-x ** 2) / (4 * (a ** 2))))
    return val


# Checking normalization function
def checknorm(func):
    norm = 0.
    for i in range(N):
        norm = (abs(func[i, :]) ** 2) + norm
    value = 1 - (dx * norm)
    return value


# Allocating full matricies
tevo1 = np.zeros((N, N), dtype=complex)
tevo2 = np.zeros((N, N), dtype=complex)
wavfunc = np.zeros((N, 1), dtype=complex)

# Allocating parts of matricies for TDMA
tevo1top = np.zeros((N-1, 1), dtype=complex)
tevo1diag = np.zeros((N, 1), dtype=complex)
tevo1bot = np.zeros((N-1, 1), dtype=complex)
# tevo2top = np.zeros((N-1, 1), dtype=complex)
# tevo2diag = np.zeros((N, 1), dtype=complex)
# tevo2bot = np.zeros((N-1, 1), dtype=complex)

# Filling tevo2
for i in range(N):
    if i == 0:
        tevo2[i, :] = [
            1 - 2j * alpha - 1j * ui if j == 0
            else 1j * alpha if j == 1
            else 0 for j in range(N)]
    if i == N:
        tevo2[i, :] = [
            1j * alpha if j == N - 1
            else 1 - 2j * alpha - 1j * ui if j == N
            else 0 for j in range(N)]
    else:
        tevo2[i, :] = [
            1j * alpha if j == i - 1 or j == i + 1
            else 1 - 2j * alpha - 1j * ui if j == i
            else 0 for j in range(N)]

# Filling the wavefunction
for i in range(N):
    wavfunc[i, :] = NormGauss(1, (i-(N/2))*dx)

# Filling in the top and bottom parts of the two matricies for the TDMA solver.
# tevo2 doesn't need to be filled for this part because tevo1 is the matrix
# that will be used in the TDMA solver.
for i in range(N-1):

    tevo1bot[i, :] = -alpha*1j
    tevo1top[i, :] = -alpha*1j
    # tevo2bot[i, :] = alpha*1j
    # tevo2top[i, :] = alpha*1j

# Filling in the middle diagonal parts of the two matricies for the TDMA solver
# Again tevo2 doesn't need to be filled this way.
for i in range(N):
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
