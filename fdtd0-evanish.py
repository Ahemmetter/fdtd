#-*- coding: utf-8 -*-

"""This program simulates the electric and magnetic field components of a Gaussian source originating at the center
of a box with an electric wall (E=0) on the left and a magnetic wall (H = 0) on the right.
The results of the FDTD simulation are displayed as an animation. Blue line: E_z, green line: H_z.

Observations: The fields propagate outwards without changing form (Gaussian) until they meet the boundaries.
At the boundaries they are reflected and move with equal but opposite velocity towards the center. On approaching the
center they overlap, the E field is canceled to nearly 0 (within numerical error) and the entire energy goes to the
H field. The animation stops at this point."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as anim

# --- Constants, Parameters and Arrays --- #

Z = 377.0                   # impedance of free space
x_length = 300              # length of the simulation area
sigma = 20.0                # width of the source
pos = 0.5 * x_length        # position of source (centered)
d = 5 * sigma               # delay
steps = int(x_length + d)   # number of steps
all_steps = np.linspace(0, x_length - 1, x_length)                              # total number of steps
E_z = np.zeros(x_length)    # empty array for E_z values
H_y = np.zeros(x_length)    # empty array for H_y values

# --- Functions --- #

def source(t, sigma, d):
    """Defines the source as a Gaussian one with width sigma."""

    return np.exp(- ((t - d) ** 2) / (2.0 * (sigma ** 2)))

def init():
    """Initializes the objects line1 and line2 as empty arrays, which will be filled with the data that should
    be plotted."""

    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,

def animate(time):
    """Main function. Calculates the field values each time step for each cell. H_y and E_z are related
    by the impedance of free space (Z). The results are written to line1 and line2 and returned."""

    H_y[:x_length - 1] = H_y[:x_length - 1] + (E_z[1:] - E_z[:x_length - 1]) / Z
    E_z[1:] = E_z[1:] + (H_y[1:] - H_y[:x_length - 1]) * Z
    E_z[pos] = E_z[pos] + source(time, sigma, d)                                # source

    line1.set_data(all_steps, E_z)                                              # values are written to line objects
    line2.set_data(all_steps, Z * H_y)
    return line1, line2,

# --- Animation --- #

fig = plt.figure()                                                              # creates figure
ax = plt.axes(xlim = (0, x_length), ylim = (-1.0, 1.0))                             # creates axes
line1, = ax.plot([], [])                                                        # plots line objects
line2, = ax.plot([], [])

anim = anim.FuncAnimation(fig, animate, init_func = init,                       # animation command
                               frames = steps, interval = 20, blit = False, repeat = False)

#anim.save('fdtd0-evanish.mp4',fps=30)                                          # optional command for saving as mp4
plt.show()