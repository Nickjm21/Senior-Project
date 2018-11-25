# File to easily share variables acorss multiple modules

# Defining variables
mass = 5.68563  # eV/(nm/fs)^2
hbar = 0.658212  # eV fs
ko = 6  # 1/nm
xo = -10  # nm
deltax = 1  # nm
deltap = hbar / (2 * deltax)  # (eV fs)/nm

L = 50  # Range of x values (nm)
N = 1000  # Number of points
dx = L / N  # Size of position steps (1/nm)
nsteps = 10  # Number of time steps
dt = 0.1  # Size of time step
alpha = (dt * hbar)/(4 * mass * (dx**2))
