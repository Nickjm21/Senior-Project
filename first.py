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

dx = 1  # Size of position steps
xrange = 10  # Range of x values
N = int(round(xrange / dx))  # Number of equations kinda
hN = int(round(xrange/2))  # Number of position steps

nsteps = 10  # Number of time steps
dt = 1  # Size of time step
pot = 1
alpha = (dt * hbar)/(4 * mass * (dx**2))
ui = (dt / (2 * hbar)) * pot  # Unsure how to implement the potential


# Function to increment in decimal values. x is initial, y is final, jump is
# the increment size.
def drange(x, y, jump):
    m = 0
    r = x
    while r < y:
        m += 1
        r = x + m * jump
        yield round(r, 5)


# Form of the initial wavefunction (Gaussian)
def Func(x):
    val = np.exp(-(x ** 2))
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

# print(tevo2)

# Filling the initial wave function (Gaussian) !This doesn't work:
for i in range(0, N):
    for k in drange(-hN, xrange * (i+1), dx):
        # print(k-hN)
        wavfunc[i, :] = Func(k - hN)

# Testing:
for i in range(N):
    wavfunc[i, :] = Func(i-hN)

print(wavfunc)
# plt.plot(abs(wavfunc))
# plt.show()

# for i in drange(-hN, hN, dx):
# print(i)

# Having a difficult time filling the initial wavefunction for non-integer
# position steps. Not allowed to assign numbers to the -4.9th element of an
# array obviously.

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

# Testing...
# print(tevo1diag)

# TDMA(tevo2bot, tevo2diag, tevo2top, wavfunc)

# np.matmul(tevo2, TDMA(tevo2bot, tevo2diag, tevo2top, wavfunc))
