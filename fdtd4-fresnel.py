import numpy as np

Z=377.0  # Free space impedance

x_length=1800  # Domain size
#Dielectric distribution
e1 = 8
e2 = 2
n1 = np.sqrt(e1)
n2 = np.sqrt(e2)

sigma = 20.0 * max(n1, n2)  # source width
pos = 0.1 * x_length  # position of source
d = sigma * 10

def source(t, sigma, d):  # definition of the source
    return np.exp(-((t - d) ** 2) / (2.0 * (sigma ** 2)))  # Gaussian distribution

def fresnel(n1, n2):
    fa = 1 - np.abs((n1 - n2) / (n1 + n2)) ** 2
    return fa

def fdtdRatio(Emax1, Emax2, Hmax1, Hmax2):
    fr = Emax2*Hmax2/(Emax1*Hmax1)
    return fr

def error(fa, fr):
    error = np.abs((fr - fa) / fa)
    return error

def fields(e1, e2):
    eps = np.ones(x_length)
    eps[:int(x_length / 2.0)] = e1
    eps[int(x_length / 2.0):] = e2

    # Left boundary
    c = 1 / np.sqrt(eps[0])
    al = (c - 1) / (c + 1)
    bl = 2 / (c + 1)
    lm, l0, lp = 0, 0, 0
    plm, pl0, plp = 0, 0, 0

    # Right boundary
    c = 1 / np.sqrt(eps[-1])
    ar = (c - 1) / (c + 1)
    br = 2 / (c + 1)
    rm, r0, rp = 0, 0, 0
    mrm, mr0, mrp = 0, 0, 0
    # Monitor points
    Emax1_x = int(2.0 * x_length / 15.0)
    Emax2_x = int(14.0 * x_length / 15.0)
    Emax1, Emax2, Hmax1, Hmax2 = 0, 0, 0, 0
    # Model
    total_steps = int(x_length * 3 + d)  # Time stepping
    # frame_interval = int(total_steps / 15.0)
    # all_steps = np.linspace(0, x_length - 1, x_length)

    # Inital field E_z and H_y is equal to zero
    E_z = np.zeros(x_length)
    H_y = np.zeros(x_length)

    for time in xrange(total_steps):

        H_y[0:x_length - 1] = H_y[0:x_length - 1] + (E_z[1:x_length] - E_z[0:x_length - 1]) / Z
        # Evaluate Mur ABC value (eq. 6.35 Taflove)
        mrp = H_y[x_length-2]
        rp = -mrm + ar * (mrp + rm) + br * (r0 + mr0)
        H_y[x_length-1] = rp
        # Cycle field values at boundary
        rm, mrm = r0, mr0
        r0, mr0 = rp, mrp

        E_z[1:x_length] = E_z[1:x_length] + Z * (H_y[1:x_length] - H_y[0:x_length - 1])/eps[0:x_length-1]
        E_z[pos] = E_z[pos] + source(time, sigma, d)
        # Evaluate Mur ABC value (eq. 6.35 Taflove)
        plp = E_z[1]
        lp = -plm + al * (plp + lm) + bl * (l0 + pl0)
        E_z[0] = lp
        # Cycle field values at boundary
        lm, plm = l0, pl0
        l0, pl0 = lp, plp

        Emax1 = max(Emax1, np.abs(E_z[Emax1_x]))
        Emax2 = max(Emax2, np.abs(E_z[Emax2_x]))
        Hmax1 = max(Hmax1, np.abs(H_y[Emax1_x]))
        Hmax2 = max(Hmax2, np.abs(H_y[Emax2_x]))

    return Emax1, Emax2, Hmax1, Hmax2


Emax1, Emax2, Hmax1, Hmax2 = fields(e1, e2)
fa = fresnel(n1, n2)
print fa
fr = fdtdRatio(Emax1, Emax2, Hmax1, Hmax2)
print fr
print error(fa, fr)
