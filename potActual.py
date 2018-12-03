import numpy as np
import config as c

omega = np.sqrt(c.ko / c.mass)


def freepart(x):
    val = 0 * x
    return val


def MM(Vo, x1, a, x):
    val = Vo * np.exp(
        -((x - x1) ** 2)
        / 2 * (a ** 2)
    )
    return val


def SHO(x):
    val = (1 / 2) * c.mass * (omega ** 2) * (x ** 2)
    return val


def box(x):
    if x < c.N/2:
        val = 0
    else:
        val = 10
    return val
