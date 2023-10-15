"""
Shoot method for finite and infinite potential well
"""
import numpy as np
from scipy.optimize import newton


def well_up(x, a, b, E):
    if a <= x <= b:
        return 0
    else:
        return E


def well_down(x, a, b, E):
    if a <= x <= b:
        return E
    else:
        return 0


def well_up_iter(x, a, b, E):
    output = []
    for element in np.nditer(x):
        if a <= element <= b:
            output.append(0)
        else:
            output.append(E)

    return np.asarray(output)


def well_down_iter(x, a, b, E):
    output = []
    for element in np.nditer(x):
        if a <= element <= b:
            output.append(E)
        else:
            output.append(0)

    return np.asarray(output)


def bvp_shoot(f, a, b, z1, z2, t, tol):
    max_iter = 30
    n = len(t)
    y1 = rku4(f, [a, z1], t)
    w1 = y1[n - 1, 0]

    for i in range(max_iter):
        y2 = rku4(f, [a, z2], t)
        w2 = y2[n - 1, 0]

        if np.abs(b - w2) < tol:
            break

        # z1 = z2
        try:
            z2 = z2 + (z2 - z1) / (w2 - w1) * (b - w2)

        except:
            print("Sumting Wong")
            break

        w1 = w2

    if np.abs(b - w2) >= tol:
        print("Maximum number of iterations exceeded. Tol fucked")

    return y2[:, 0]


def symp_pefrl(f, psi0, x, V, E):
    """
    Solves 2. order IVP using "Position Extended Forest-Ruth Like" algorithm of Omelyan et al.

    INPUTS:
    f: function to equal dpsi^2/dx^2
    x0: initial value for f(t[0])
    v0: initial value for df/dx(v[0])
    x: array of positions to evaluate function

    OUTPUT:
    output: 2D Array of evaluated psi(x) values and evaluated df/dx(x, t) values
    """
    xi = 0.1786178958448091
    lam = -0.2123418310626054
    chi = -0.6626458266981849e-1
    # h = t[1] - t[0]
    x_prev = x[0]
    x0, v0 = psi0
    x_pre = x0
    v = v0
    output = []
    c = 0  # Debug index counter for V to prevent 40k dim array
    for pos in x:

        h = pos - x_prev
        x_pre += xi * h * v
        v += (1 - 2*lam)*0.5*h*f(x_pre, pos, V[c], E)
        x_pre += chi * h * v
        v += lam*h*f(x_pre, pos, V[c], E)
        x_pre += (1 - 2 * (chi + xi)) * h * v
        v += lam*h*f(x_pre, pos, V[c], E)
        x_pre += chi * h * v
        v += (1 - 2*lam)*0.5*h*f(x_pre, pos, V[c], E)
        x_pre += xi * h * v
        x_prev = pos
        c += 1
        output.append([x_pre, v])

    return np.column_stack(np.array(output))


def schrodinger(y, r, V, E):
    psi, phi = y
    dphidx = [phi, (V - E) * psi]
    return np.asarray(dphidx)


def schrodinger2(y, r, V, E):
    return (V - E) * y


def rku4(f, x0, t, V, E):
    n = len(t)
    x = np.array([x0] *n)
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = h * f(x[i], t[i], V[i], E)
        k2 = h * f(x[i] + 0.5 * k1, t[i] + 0.5 * h, V[i], E)
        k3 = h * f(x[i] + 0.5 * k2, t[i] + 0.5 * h, V[i], E)
        k4 = h * f(x[i] + k3, t[i+1], V[i], E)
        x[i+1] = x[i] + (k1 + 2.0 * (k2 + k3) + k4) / 6.0

    return x


def rk4_psi(f, psi0, x, V, E):
    """RK4 method that also takes V and E parameters. Made to solve IVP phi'=f(psi, x) with psi(x[0])=psi0."""
    n = len(x)
    psi = np.array([psi0]*n)
    for i in range(n - 1):
        h = x[i+1] - x[i]
        k1 = h * f(psi[i], x[i], V[i], E)
        k2 = h * f(psi[i] + 0.5*k1, x[i] + 0.5*h, V[i], E)
        k3 = h * f(psi[i] + 0.5*k2, x[i] + 0.5*h, V[i], E)
        k4 = h * f(psi[i] + k3, x[i+1], V[i], E)
        psi[i+1] = psi[i] + (k1 + 2*(k2 + k3) + k4)/6

    return psi


def shoot_psi(f, psi0, x, V, E_arr):
    """"Shooting method: find zeroes of Schrödinger equation f with potential V for energies in array E_arr"""
    psi_right = []
    for energy in E_arr:
        psi = symp_pefrl(f, psi0, x, V, energy)
        dim = np.shape(psi)[1]
        psi_right.append(psi[0][dim - 1])

    return np.array(psi_right)


def one_shot(E, f, psi0, x, V):  # Reordered inputs for scipy newton
    """Same as shoot_psi but only for one value of energy"""
    psi = symp_pefrl(f, psi0, x, V, E)
    return psi[0][np.shape(psi)[1] - 1]


def find_zeroes(rightbound_vals):
    """Amazing method found online. Find zero crossing due to sign change in input array."""
    return np.where(np.diff(np.signbit(rightbound_vals)))[0]


def optimize_energy(f, psi0, x, V, E_arr):
    shoot_try = shoot_psi(f, psi0, x, V, E_arr)  # Should now be a 2D array with psi and psi' right bound values
    crossings = find_zeroes(shoot_try)
    energy_list = []
    for cross in crossings:
        # Use Newton-Raphson method to find zero of function
        energy_list.append(newton(one_shot, E_arr[cross], args=(f, psi0, x, V)))

    return np.array(energy_list)


def normalize(wavefunction, pos, width):
    """Some sort of really rough normalization"""
    return wavefunction/np.max(wavefunction[pos:(pos + width)])


def normalize_all(wavefunction):
    """Some sort of really rough normalization"""
    return wavefunction/np.max(wavefunction)


def normalize_range(wavefunction, start, end):
    """Some sort of really rough normalization"""
    dim = len(wavefunction)
    a = int(start*dim)
    b = int(end*dim)
    return wavefunction/np.max(wavefunction[a:b])


def potwell(psi_init, upper, h_):
    """Solves infinite potential well numerically and analytically. Also returns eigenenergies."""
    x_arr_ipw = np.arange(0, 1+h_, h_)
    V_ipw = np.zeros(len(x_arr_ipw))
    E_arr = np.arange(1, upper, 5)  # Set initial guesses for eigen energies
    eigE = optimize_energy(schrodinger2, psi_init, x_arr_ipw, V_ipw, E_arr)
    ipw_output = []
    for energy in eigE:
        # Get numerical solution for each eigen energy
        out = symp_pefrl(schrodinger2, psi_init, x_arr_ipw, V_ipw, energy)
        ipw_output.append(normalize_all(out[0, :]))

    # Get analytical solutions
    k = np.arange(1, len(eigE) + 1)
    ipw_solve_analytical = []
    for kk in k:
        ipw_solve_analytical.append(np.sin(kk*np.pi*x_arr_ipw))

    return x_arr_ipw, np.array(ipw_output), np.array(ipw_solve_analytical), eigE


def one_shot_fin(E, f, psi0, x, V):  # Reordered inputs for scipy newton
    """Same as shoot_psi but only for one value of energy"""
    psi = symp_pefrl(f, psi0, x, V, E)
    dim = np.shape(psi)[1]
    return psi[0, dim - 1]


def shoot_psi_fin(f, psi0, x, V, E_arr):
    """"Shooting method: find zeroes of Schrödinger equation f with potential V for energies in array E_arr"""
    psi_right = []
    mult_psi = []
    for energy in E_arr:
        psi = symp_pefrl(f, psi0, x, V, energy)
        dim = np.shape(psi)[1]
        psi_right.append(psi[0, dim - 1])
        mult_psi.append(psi)

    return np.array(psi_right)


def optimize_energy_fin(f, psi0, x, V, E_arr):
    shoot_try = shoot_psi_fin(f, psi0, x, V, E_arr)
    # TODO: Fix boundry value here for 0
    crossings = find_zeroes(shoot_try)
    energy_list = []
    for cross in crossings:
        # Use Newton-Raphson method to find zero of function
        energy_list.append(newton(one_shot_fin, E_arr[cross], args=(f, psi0, x, V)))

    return np.array(energy_list)


def finpotwell(psi_init, upper, depth, h_):
    """Solves finite potential well numerically. Also returns eigenenergies."""
    x_arr_fpw = np.arange(-10, 10 + h_, h_)
    dim = len(x_arr_fpw)
    pos = int(dim//2.2)
    width = int(2*(dim/2 - pos))
    V_fpw = np.zeros(dim)
    V_fpw[:pos] = depth
    V_fpw[(pos + width):] = depth
    E_arr = np.arange(1, upper, 5)  # Set initial guesses for eigen energies
    eigE = optimize_energy_fin(schrodinger2, psi_init, x_arr_fpw, V_fpw, E_arr)
    fpw_output = []
    for energy in eigE:
        # Get numerical solution for each eigen energy
        out = symp_pefrl(schrodinger2, psi_init, x_arr_fpw, V_fpw, energy)
        fpw_output.append(normalize(out[0, :], pos, width))

    return x_arr_fpw, np.array(fpw_output), eigE, V_fpw


# TODO: Fuck it, try LHO
def lho(psi_init, upper, h_):
    x_arr_lho = np.arange(-10, 10 + h_, h_)
    V_lho = x_arr_lho**2
    E_arr = np.arange(1, upper, 5)  # Set initial guesses for eigen energies
    eigE = optimize_energy_fin(schrodinger2, psi_init, x_arr_lho, V_lho, E_arr)
    fpw_output = []
    for energy in eigE:
        # Get numerical solution for each eigen energy
        out = symp_pefrl(schrodinger2, psi_init, x_arr_lho, V_lho, energy)
        fpw_output.append(normalize_range(out[0, :], 1/4, 3/4))

    return x_arr_lho, np.array(fpw_output), eigE, V_lho

