import numpy as np
import matplotlib.pyplot as plt
from src_shoot import *
import scipy.integrate

# Params
V = -1000
E = -200


def schrod(x, t):
    psi, phi = x
    dphidx = [phi, (V - E) * psi]
    return np.asarray(dphidx)


def symp_schrod(x, t, V, E):
    return (V - E) * x


psi_0 = 0
phi_0 = 1  # psi' = phi
psi_init = np.array([psi_0, phi_0])
h_ = 1/2000  # Step size
a, b = (0.0, 1)
x0 = 1.0
n = 2000
x = np.linspace(a, b, n)
V = np.zeros(len(x))
# sci = np.column_stack(scipy.integrate.odeint(schrod, psi_init, x))
sol = symp_pefrl(symp_schrod, (0, 1), x, V, 100)
# sci = np.array([float(value) for value in np.nditer(sci_pre)])
# plt.plot(x, sci[1])
plt.plot(x, sol[0])
plt.show()
