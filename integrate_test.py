import numpy as np
import matplotlib.pyplot as plt
from src_shoot import *
import scipy.integrate


def f(x, t, V, E):
    return t*np.sin(t)


def f2(x, t):
    return t*np.sin(t)


a, b = (0.0, 20.0)
x0 = 1.0

n = 2000
V = np.zeros(n)
t = np.linspace(a, b, n)
exact = np.sin(t) - t*np.cos(t) + x0
sci_py = scipy.integrate.odeint(f2, x0, t)
solve = rk4_psi(f, x0, t, V, 0)


plt.plot(t, solve, label="rk4_psi")
plt.plot(t, exact, label="Exact")
plt.legend()
plt.show()

plt.plot(t, np.abs(solve - exact))
plt.show()