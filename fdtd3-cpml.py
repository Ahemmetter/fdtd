#-*- coding: utf-8 -*-

"""This program simulates the behaviour of a CPML (convolution perfectly matched layer).

Observation: The wave seems to be absorbed by the walls."""

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
steps = int(round((x_length + d) * n))                    # number of steps
all_steps = np.linspace(0, x_length-1, x_length)    # list of all steps
E_z = np.zeros(x_length)                            # creates empty list of E and H values
H_y = np.zeros(x_length)

# for the CPML layer:

dx = 1.0
R0 = 1e-5
m = 2.85  # Order of polynomial grading
pml_width = 5.0
sxmax = -(m+1)*np.log(R0)/2/Z/(pml_width*dx)
sx = np.zeros(x_length)
sxm = np.zeros(x_length)
Phx = np.zeros(x_length)
Pex = np.zeros(x_length)

for mm in xrange(int(pml_width)):

    sx[mm+1] = sxmax*((pml_width-mm-0.5)/pml_width)**m
    sxm[mm] = sxmax*((pml_width-mm)/pml_width)**m  # Shifted to the right
    sx[x_length-mm-1] = sxmax*((pml_width-mm-0.5)/pml_width)**m
    sxm[x_length-mm-1] = sxmax*((pml_width-mm)/pml_width)**m

aex = np.exp(-sx*Z)-1
bex = np.exp(-sx*Z)
ahx = np.exp(-sxm*Z)-1
bhx = np.exp(-sxm*Z)

# --- Functions --- #

def source(t, sigma, d):                       # definition of the source
    return np.exp(-((t-d)**2)/(2.0*(sigma**2))) # Gaussian distribution

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return line1, line2,

def animate(time):

    Phx[:x_length - 1] = bhx[:x_length - 1] * Phx[:x_length - 1] + ahx[:x_length - 1] * (E_z[1:x_length] - E_z[:x_length - 1])
    H_y[:x_length - 1] = H_y[:x_length - 1] + (E_z[1:x_length] - E_z[:x_length - 1]) / Z + Phx[:x_length - 1] / Z
    Pex[1:x_length] = bex[1:x_length] * Pex[1:x_length] + aex[1:x_length] * (H_y[1:x_length] - H_y[:x_length - 1])
    E_z[1:x_length] = E_z[1:x_length] + (H_y[1:x_length] - H_y[:x_length - 1]) * Z / eps[:x_length - 1] + Pex[1:x_length] * Z / eps[:x_length - 1]
    E_z[pos] = E_z[pos] + source(time, sigma, d)

    line1.set_data(all_steps, E_z)
    line2.set_data(all_steps, Z * H_y)
    return line1, line2,

# --- Animation --- #

fig = plt.figure()  # creates figure
ax = plt.axes(xlim = (0, x_length), ylim = (-epsilon, epsilon)) # creates axes
line1, = ax.plot([], [])                            # creates empty line objects
line2, = ax.plot([], [])

anim = anim.FuncAnimation(fig, animate, init_func = init,
                               frames = steps, interval = 20, blit = False, repeat = False)

#anim.save('fdtd3-cpml.mp4',fps=30)
plt.show()