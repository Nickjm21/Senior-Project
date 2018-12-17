# File to easily share variables acorss multiple modules
import numpy as np

# Defining variables
mass = 5.68563  # eV/(nm/fs)^2
hbar = 0.658212  # eV fs
ko = 6  # 1/nm
xo = -10  # nm
deltax = 3  # nm
deltap = hbar / (2 * deltax)  # (eV fs)/nm

xmin = -50.  # (nm)
xmax = 50.  # (nm)
N = 2500  # Number of points
dx = (xmax - xmin) / N  # Size of position steps (1/nm)
nsteps = 1  # Number of time steps
dt = 50 / 300  # Size of time step
alpha = (dt * hbar)/(4 * mass * (dx**2))

Vo = 1.6  # eV
a = 1.5  # nm
x1 = -5  # nm
potmin = 0  # nm
potmax = 2  # nm


def deltaE(k):
    val = (
        np.sqrt(
        ((hbar ** 2) * (8 * (k ** 2) * (deltax ** 2) + 1))
        / ((mass ** 2) * (deltax ** 4))) / (4 * np.sqrt(2)))
    return val

deltaE(10)
