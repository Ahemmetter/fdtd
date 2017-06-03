#-*- coding: utf-8 -*-

"""This program simulates the behaviour of a Mur absorbing boundary condition.

Observation: The wave seems to travel outside of the simulated area."""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as anim

# --- Constants, Parameters and Arrays --- #

Z = 377.0                   # impedance of free space
x_length = 300              # length of the simulation area
epsilon = 1                 # dielectric permittivity
eps = np.zeros(x_length)    # empty array for permittivity values
eps[:] = epsilon            # set permittivity for each cell
n = np.sqrt(epsilon)        # refractive index
sigma = 20.0 * n            # source width
pos = 0.5 * x_length        # position of source
d = 5 * sigma               # delay
c = 1.0 / n                 # speed of light in the medium
a = (c-1)/(c+1)             # ABC constants
b = 2/(c + 1)
steps = int(round(x_length + d))                    # number of steps
all_steps = np.linspace(0, x_length-1, x_length)    # list of all steps
E_z = np.zeros(x_length)                            # creates empty list of E and H values
H_y = np.zeros(x_length)

# --- Functions --- #

def source(t, sigma, d):
    """Defines the source as a Gaussian one."""

    return np.exp( - ((t - d) ** 2) / (2.0 * (sigma ** 2)))

def init():
    """Initializes line objects."""

    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,

def animate(time, lm=0, l0=0, lp=0, plm=0, pl0=0, plp=0, rm=0, r0=0, rp=0, mrm=0, mr0=0, mrp=0):
    """Main function. Calculates the field values for each time step and applies boundary conditions."""

    # H field:

    H_y[:x_length - 1] = H_y[:x_length - 1] + (E_z[1:x_length] - E_z[:x_length - 1]) / Z
    mrp = H_y[-2]
    rp = - mrm + a * (mrp + rm) + b * (r0 + mr0)
    H_y[-1] = rp
    rm = r0
    mrm = mr0
    r0 = rp
    mr0 = mrp

    # E field:

    E_z[1:x_length] = E_z[1:x_length] + (H_y[1:x_length] - H_y[:x_length - 1]) * Z
    E_z[pos] = E_z[pos] + source(time, sigma, d)
    plp = E_z[1]
    lp = - plm + a * (plp + lm) + b * (l0 + pl0)
    E_z[0] = lp
    lm = l0
    plm = pl0
    l0 = lp
    pl0 = plp

    line1.set_data(all_steps, E_z)
    line2.set_data(all_steps, Z * H_y)
    return line1, line2,

# --- Animation --- #

fig = plt.figure()                                      # creates figure
ax = plt.axes(xlim = (0, x_length), ylim = (-1.0, 1.0)) # creates axes
line1, = ax.plot([], [])                                # creates empty line objects
line2, = ax.plot([], [])

anim = anim.FuncAnimation(fig, animate, init_func = init,
                               frames = steps, interval = 20, blit = False, repeat = False)

#anim.save('fdtd2-mur.mp4',fps=30)                       # optional command to save mp4
plt.show()