import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as anim

# --- Constants --- #

Z = 377.0                           # impedance of free space
x_length = 100                      # domain size
l = int(x_length / 5.0)             # wavelength
factor = 2.0                        # for slab lambda/2
e1 = 1                              # permittivity of surrounding material
e2 = 2                              # permittivity of slab
n1 = np.sqrt(e1)                    # refractive indices
n2 = np.sqrt(e2)

# --- Epsilon --- #

print("Epsilon in slab was set to be %f"%e2)
n2 = (l/factor)/int(l/n2/factor)    # corrected refractive index
e2 = n2 ** 2                        # corrected e2
thickness = l / n2 / factor         # slab thickness
print("Corrected epsilon in slab is %f (to fit exactly to space step)"%e2)
er = np.zeros(x_length)
e = np.zeros(x_length)
er[:] = e1
e[:] = e1
e[int((x_length / 2.0) - thickness / 2.0):int(x_length / 2.0 + thickness / 2.0)] = e2

c = 1/np.sqrt(e[0])
cref = c

sigma = l / 3.0  # source width
pos = 0.1 * x_length  # position of source
d = sigma * 5

def source(t, sigma, d):  # definition of the source
    return np.exp(-((t - d) ** 2) / (2.0 * (sigma ** 2)))
    # for continuous source
    # if t > d:
    #     amp = 1.0
    # amp = 1.0
    # return amp/ n1 * np.sin(2 * np.pi * t * cref * n1 / l)

steps = int(x_length*4.0+d)  # Time stepping
frame_interval = int(steps/35.0)

E_z = np.zeros(x_length)
H_y = np.zeros(x_length)
refez = np.zeros(x_length)
refhy = np.zeros(x_length)
x = np.arange(0, x_length, 1)

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    line6.set_data([], [])
    return line1, line2, line3, line4, line5, line6,

def animate(time):

    # H field
    H_y[-1] = H_y[-2]
    H_y[:-1] = H_y[:-1] + (E_z[1:] - E_z[:-1]) / Z

    # H field reference
    refhy[-1] = refhy[-2]
    refhy[:-1] = refhy[:-1] + (refez[1:] - refez[:-1]) / Z

    # E field
    E_z[0] = E_z[1]
    E_z[1:] = E_z[1:] + (H_y[1:] - H_y[:-1]) * Z / e[1:]
    E_z[pos] += source(time, sigma, d)

    # E field reference
    refez[0] = refez[1]
    refez[1:] = refez[1:] + (refhy[1:] - refhy[:-1]) * Z / er[1:]
    refez[pos] += source(time, sigma, d)

    line1.set_data(x, E_z)
    line2.set_data(x, Z * H_y)
    line3.set_data(x, refez)
    line4.set_data(x, Z * refhy)
    line5.set_data(x, abs(E_z - refez))
    line6.set_data(x, Z * abs(H_y - refhy))
    return line1, line2, line3, line4, line5, line6,

# --- Animation --- #

fig, (a0, a1) = plt.subplots(2, sharex = True, gridspec_kw = {"height_ratios": [3, 1]})
a0.set_xlim([0, x_length])
a0.set_ylim([-1.0, 1.0])
a1.set_ylim([0, 1.0])
line1, = a0.plot([], [])                # creates empty line objects
line2, = a0.plot([], [])
line3, = a0.plot([], [])
line4, = a0.plot([], [])
line5, = a1.plot([], [])
line6, = a1.plot([], [])
a0.axvspan(x_length/2.0 - thickness, x_length/2.0 + thickness, alpha=0.3, color='green')
a1.axvspan(x_length/2.0 - thickness, x_length/2.0 + thickness, alpha=0.3, color='green')

anim = anim.FuncAnimation(fig, animate, init_func=init,
                          frames=steps, interval=30, blit=False, repeat=False)

# anim.save('fdtd5-slab.mp4',fps=30)
plt.show()