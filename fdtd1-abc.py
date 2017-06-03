#-*- coding: utf-8 -*-

"""This program simulates a simple absorbing boundary condition (ABC) by forcing the field values immediately
in front of the boundaries to the value preceding it.

Observation: The wave seems to travel outside of the simulated area."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as anim

# --- Constants, Parameters and Arrays --- #

Z = 377.0                   # impedance of free space
x_length = 300              # length of the simulation area
epsilon = 1                 # dielectric permittivity
eps = np.zeros(x_length)    # empty array for permittivity values
eps[:] = epsilon            # set the permittivity for each cell
n = np.sqrt(epsilon)        # refractive index for non-magnetic material
sigma = 20.0 * n            # source width
pos = 0.5 * x_length        # position of source (centered)
d = 5 * sigma               # delay
E_z = np.zeros(x_length)    # empty list for E_z values
H_y = np.zeros(x_length)    # empty list for H_y values
steps = int(round(x_length * n))                                                # number of steps, rounded for odd n
all_steps = np.linspace(0, x_length - 1, x_length)                              # list of all steps

# --- Functions --- #

def source(t, sigma, d):
    """Defines the source as a Gaussian one with width sigma."""

    return np.exp(- ((t - d) ** 2) / (2.0 * (sigma ** 2)))

def init():
    """ Initializes the line objects to be filled with values."""

    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,

def animate(time):
    """Main function. Simulates the propagation and absorption of the EM wave packet
    and writes it to the line objects."""

    H_y[x_length - 1] = H_y[x_length - 2]
    H_y[:x_length - 1] = H_y[:x_length - 1] + (E_z[1:x_length] - E_z[:x_length - 1]) / Z
    E_z[0] = E_z[1]
    E_z[1:x_length] = E_z[1:x_length] + (H_y[1:x_length] - H_y[0:x_length - 1]) * Z
    E_z[pos] = E_z[pos] + source(time, sigma, d)

    line1.set_data(all_steps, E_z)
    line2.set_data(all_steps, Z * H_y)
    return line1, line2,

# --- Animation --- #

fig = plt.figure()                                                              # creates figure
ax = plt.axes(xlim = (0, x_length), ylim = (-1.0, 1.0))                         # creates axes
line1, = ax.plot([], [])                                                        # creates empty line objects
line2, = ax.plot([], [])

anim = anim.FuncAnimation(fig, animate, init_func = init,                       # animation command
                               frames = steps, interval = 20, blit = False, repeat = False)

#anim.save("fdtd1-abc.mp4",fps=30)                                               # optional command to save mp4
plt.show()