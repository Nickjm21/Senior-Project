import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time_evolve as tE

fig, ax = plt.subplots()

x = np.arange(1, 500, 1)
line, = ax.plot(x, tE.timeEvo(x, tE.wavfunc))

line,


def func(frame, num, wavefunction):
    val = tE.timeEvo(num, wavefunction)
    return val


def init():  # only required for blitting to give a clean slate.
    line.set_ydata([np.nan] * len(x))
    return line,


def animate(i):
    line.set_data(tE.timeEvo(i, tE.wavfunc))  # update the data.
    return line,


ani = animation.FuncAnimation(
    fig, animate, init_func=init, interval=2, blit=True, save_count=50)


animation.FuncAnimation(
    fig, func(i, tE.wavfunc), frames=None, init_func=None, fargs=None, save_count=None)
# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()
