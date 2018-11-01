import numpy
import numpy as np
import matplotlib
import matrix
import matplotlib.pyplot as plt
import decimal

TDMA = matrix.TDMAsolver
# Defining variables
mass = 1
hbar = 1

dx = 0.1  # Size of position steps
range = 10  # Range of x values
N = int(round(range / dx))  # Number of equations kinda
hN = int(round(range/2))  # Number of position steps

nsteps = 10  # Number of time steps
dt = 1  # Size of time step
pot = 1

alpha = (dt * hbar)/(4 * mass * (dx**2))
ui = (dt / (2 * hbar)) * pot  # Unsure how to implement the potential


# Increments in decimal values
def drange(x, y, jump):
    m = 0
    r = x
    while r < y:
        m += 1
        r = x + m * jump
        yield round(r, 5)


# Testing


# Form of the initial wavefunction
def Func(x):
    val = np.exp(- (x ** 2))
    return(val)


# Allocating full matricies
tevo1 = np.zeros((N, N), dtype=complex)
tevo2 = np.zeros((N, N), dtype=complex)
wavfunc = np.zeros((N, 1), dtype=complex)

# Allocating parts of matricies for TDMA

tevo1top = np.zeros((N-1, 1), dtype=complex)
tevo1diag = np.zeros((N, 1), dtype=complex)
tevo1bot = np.zeros((N-1, 1), dtype=complex)
tevo2top = np.zeros((N-1, 1), dtype=complex)
tevo2diag = np.zeros((N, 1), dtype=complex)
tevo2bot = np.zeros((N-1, 1), dtype=complex)

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

# Filling remainder of matricies
for j in range(0, N):
    for k in drange(j, range * (j+1), dx):
        wavfunc[j, :] = Func(k - hN)

print(wavfunc)

plt.plot(abs(wavfunc))
plt.show()


for i in range(N-1):

    tevo1bot[i, :] = -alpha*1j
    tevo1top[i, :] = -alpha*1j
    tevo2bot[i, :] = alpha*1j
    tevo2top[i, :] = alpha*1j

for i in range(N):

    tevo1diag[i, :] = 1 + 2j*alpha + 1j*ui
    tevo2diag[i, :] = 1 - 2j*alpha - 1j*ui

print(tevo1diag)

TDMA(tevo2bot, tevo2diag, tevo2top, wavfunc)

np.matmul(tevo2, TDMA(tevo2bot, tevo2diag, tevo2top, wavfunc))

plt.xlim(0, 8)
plt.plot(abs(wavfunc))
plt.show()


for a in range(10):
    print(a)
