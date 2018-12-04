# File to easily share variables acorss multiple modules

# Defining variables
mass = 5.68563  # eV/(nm/fs)^2
hbar = 0.658212  # eV fs
ko = 6  # 1/nm
xo = -10  # nm
deltax = 1  # nm
deltap = hbar / (2 * deltax)  # (eV fs)/nm

xmin = -25.  # (nm)
xmax = 25.  # (nm)
N = 100  # Number of points
dx = (xmax - xmin) / N  # Size of position steps (1/nm)
nsteps = 1  # Number of time steps
dt = 38 / 1  # Size of time step
alpha = (dt * hbar)/(4 * mass * (dx**2))

Vo = 1.2  # eV
a = 2  # nm
x1 = -5  # nm
