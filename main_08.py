import numpy as np
import cmasher as cmr
import matplotlib.pyplot as plt
import matplotlib.colors
from src_shoot import *
from src_fd import *
from bvp import *
from scipy.linalg import solve_banded
from matplotlib.animation import ArtistAnimation, FuncAnimation


# Schrodinger parameters
# E = 5
# V = -10
# a = -5
# b = 5
# k = 1
#
#
# def schrodinger_up(psi, x):
#     """
#     Schrödingers equation for well potential
#     """
#     k = 1  # hbar^2/2m
#     return np.array([psi[1], 1/k * (well_up(x, a, b, V) - E) * psi[0]])
#
#
# def schrodinger_down(psi, x):
#     """
#     Schrödingers equation for well potential
#     """
#     k = 1  # hbar^2/2m
#     return np.array([psi[1], 1/k * (well_down(x, a, b, V) - E) * psi[0]])


def error_with_dim(start, stop, state):
    """Measure average error with larger and larger matrix sizes"""
    n_range = np.arange(start, stop)
    output = []
    x_range = []
    img = []
    for n in n_range:
        x, eigE, eigPsi = fd_ipw(n, 1, 0)
        exact = analytic(x, state + 1)
        output.append(np.average(get_err(eigPsi[:, state], exact)))
        #output.append(eigPsi[:, state])
        x_range.append(x)
        #if np.average(eigPsi[:, state]) > 0:  # Fix some orientations being down
        #    img.append([plt.plot(x, eigPsi[:, state], c="#0091AD"), ttl])
        #else:
        #    img.append([plt.plot(x, -eigPsi[:, state], c="#0091AD"), ttl])

    return n_range, np.array(x_range), np.array(output), np.array(img)


def test(x, t):
    return np.array([x[1], x[0] + 4*np.exp(t)])


def test_exact(t):
    return np.exp(t) * (1 + 2*t)


# Test plots
# a = 0
# b = 0.5
# tol = 10**-6
# z1 = np.exp(a)
# z2 = 2*np.exp(b)
# t_space = np.linspace(a, b, 300)
# sol = shoot(test, z1, z2, z1, z2, t_space, tol)
# sol2 = fd(4*np.exp(t_space), 1, 0, t_space, z1, z2)
# exact = test_exact(t_space)
#
# fig, (ax1, ax2) = plt.subplots(1, 2)
# ax1.plot(t_space, sol, label="Shoot", c="#B7094C")
# ax1.plot(t_space, sol2, label="FD", c="#723C70")
# ax1.plot(t_space, exact, label="Točna rešitev", c="#0091AD")
# ax1.legend()
# ax2.plot(t_space, np.abs(sol - exact), label="Shoot", c="#B7094C")
# ax2.plot(t_space, np.abs(sol2 - exact), label="FD", c="#723C70")
# ax2.legend()
#
# plt.suptitle(r"Poskusno reševanje $x'' = x + \exp{t}$")
# plt.yscale("log")
# ax1.set_xlabel("x")
# ax1.set_ylabel(r"$f(x)$")
# ax2.set_xlabel("x")
# ax2.set_ylabel(r"$|f(x) - exact|$")
# ax2.yaxis.tick_right()
# plt.show()

# Set initial conditions and parameters
psi_0 = 0
phi_0 = 1  # psi' = phi
psi_init = np.array([psi_0, phi_0])
h_ = 1/2000  # Step size
upper = 500
depth = 100
t = np.arange(-10, 10 +h_, h_)

# ipw_x, ipw_num, ipw_ana, ipw_E = potwell(psi_init, upper, h_)
# n_solutions = np.shape(ipw_ana)[0]
# print("Found {} solutions!".format(n_solutions))
# print(ipw_E)

# # Plot
# plt.plot(ipw_x, ipw_num[0, :], label=r"$\psi_1$, $E_1 =$ {}".format(round(ipw_E[0], 2)), c="#B7094C")
# plt.plot(ipw_x, ipw_num[1, :], label=r"$\psi_2$, $E_2 =$ {}".format(round(ipw_E[1], 2)), c="#892B64")
# plt.plot(ipw_x, ipw_num[2, :], label=r"$\psi_3$, $E_3 =$ {}".format(round(ipw_E[2], 2)), c="#723C70")
# plt.plot(ipw_x, ipw_num[3, :], label=r"$\psi_4$, $E_4 =$ {}".format(round(ipw_E[3], 2)), c="#5C4D7D")
# plt.plot(ipw_x, ipw_num[4, :], label=r"$\psi_5$, $E_5 =$ {}".format(round(ipw_E[4], 2)), c="#2E6F95")
# plt.plot(ipw_x, ipw_num[5, :], label=r"$\psi_6$, $E_6 =$ {}".format(round(ipw_E[5], 2)), c="#0091AD")
# plt.title("Lastne funkcije neskončne potencialne jame")
# plt.xlabel("x/L")
# plt.ylabel(r"$\psi_n(x)$")
# plt.legend()
# plt.show()
#
# # Abs error plot
# plt.plot(ipw_x, get_err(ipw_num[0, :], ipw_ana[0, :]), label=r"$\psi_1$", c="#B7094C")
# plt.plot(ipw_x, get_err(ipw_num[1, :], ipw_ana[1, :]), label=r"$\psi_2$", c="#892B64")
# plt.plot(ipw_x, get_err(ipw_num[2, :], ipw_ana[2, :]), label=r"$\psi_3$", c="#723C70")
# plt.plot(ipw_x, get_err(ipw_num[3, :], ipw_ana[3, :]), label=r"$\psi_4$", c="#5C4D7D")
# plt.plot(ipw_x, get_err(ipw_num[4, :], ipw_ana[4, :]), label=r"$\psi_5$", c="#2E6F95")
# plt.plot(ipw_x, get_err(ipw_num[5, :], ipw_ana[5, :]), label=r"$\psi_6$", c="#0091AD")

# plt.title("Absolutna napaka numeričnih rešitev")

# plt.yscale("log")
# plt.ylabel(r"$|num - analytical|$")
# plt.xlabel("x/L")
# plt.legend(loc="upper right")
# plt.show()

# fpw_x, fpw_num, fpw_E, fpw_V = finpotwell(psi_init, upper, depth, h_)
# n_solutions = np.shape(fpw_E)[0]
# print("Found {} solutions!".format(n_solutions))
# # fpw_num = np.column_stack(fpw_num)
# plt.plot(fpw_x, 50*fpw_num[0], label=r"$\psi_1$, $E_1=$" + str(round(fpw_E[0], 2)), c="#B7094C")   # c="#00BFB2")
# plt.plot(fpw_x, 50*fpw_num[1], label=r"$\psi_2$, $E_2=$" + str(round(fpw_E[1], 2)), c="#892B64")   # c="#1A5E63")
# plt.plot(fpw_x, 50*fpw_num[2], label=r"$\psi_3$, $E_3=$" + str(round(fpw_E[2], 2)), c="#723C70")   # c="#028090")
# plt.plot(fpw_x, 50*fpw_num[3], label=r"$\psi_4$, $E_4=$" + str(round(fpw_E[3], 2)), c="#5C4D7D")   # c="#EFBDEB")
# plt.plot(fpw_x, 50*fpw_num[4], label=r"$\psi_5$, $E_5=$" + str(round(fpw_E[4], 2)), c="#2E6F95")
# plt.plot(fpw_x, 50*fpw_num[5], label=r"$\psi_6$, $E_6=$" + str(round(fpw_E[5], 2)), c="#0091AD")
#
# plt.plot(fpw_x, fpw_V, c="#151729")
#
# plt.title("Vezana stanja končne potencialne jame z globino 100 a.u.")
# plt.xlabel("x")
# plt.ylabel(r"$\psi_n(x)$")
# plt.ylim(-100, 125)
# plt.xlim(-2, 2)
# plt.legend()
# plt.show()

# FD method IPW plot
# dim = 1000
# x, eigE, eigPsi = fd_ipw(dim, 1, 0)
# plt.plot(x, eigPsi[:, 0], label=r"$\psi_1$, $E_1=$" + str(round(eigE[0], 2)), c="#B7094C")
# plt.plot(x, eigPsi[:, 1], label=r"$\psi_2$, $E_2=$" + str(round(eigE[1], 2)), c="#892B64")
# plt.plot(x, eigPsi[:, 2], label=r"$\psi_3$, $E_3=$" + str(round(eigE[2], 2)), c="#723C70")
# plt.plot(x, eigPsi[:, 3], label=r"$\psi_4$, $E_4=$" + str(round(eigE[3], 2)), c="#5C4D7D")
# plt.plot(x, eigPsi[:, 4], label=r"$\psi_5$, $E_5=$" + str(round(eigE[4], 2)), c="#2E6F95")
# plt.plot(x, eigPsi[:, 5], label=r"$\psi_6$, $E_6=$" + str(round(eigE[5], 2)), c="#0091AD")
#
# plt.title("Vezana stanja neskočne potencialne jame z FD")
# plt.xlabel("x")
# plt.ylabel(r"$\psi_n(x)$")
# plt.legend(loc="upper right")
# plt.show()

# FD method IPW abs error plot
# plt.plot(x, get_err(eigPsi[:, 0], analytic(x, 1)), label=r"$\psi_1$ abs. err.", c="#B7094C")
# plt.plot(x, get_err(eigPsi[:, 1], analytic(x, 2)), label=r"$\psi_2$ abs. err.", c="#892B64")
# plt.plot(x, get_err(eigPsi[:, 2], analytic(x, 3)), label=r"$\psi_3$ abs. err.", c="#723C70")
# plt.plot(x, get_err(eigPsi[:, 3], analytic(x, 4)), label=r"$\psi_4$ abs. err.", c="#5C4D7D")
# plt.plot(x, get_err(eigPsi[:, 4], analytic(x, 5)), label=r"$\psi_5$ abs. err.", c="#2E6F95")
# plt.plot(x, get_err(eigPsi[:, 5], analytic(x, 6)), label=r"$\psi_6$ abs. err.", c="#0091AD")
# plt.title("Absolutne napake rešitev z FD in n = {}".format(dim))
# plt.xlabel("x")
# plt.ylabel(r"$|num - analytic|$")
# plt.yscale("log")
# plt.legend()
# plt.show()

# FD method IPW abs error vs. matrix dim plot
# print(np.average(np.abs(-eigPsi[:, 0] - analytic(x, 1))))
# plt.plot(x, -eigPsi[:, 0])
# plt.plot(x, analytic(x, 1))
# plt.show()

# FD method FPW plot
# dim = 10000
# x_fpw, eigE_fpw, eigPsi_fpw, V_fpw = fd_fpw(dim, 10, 100)
# n_solutions = np.shape(eigE_fpw)[0]
# print("Found {} solutions!".format(n_solutions))
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 0], label=r"$\psi_1$, $E_1=$" + str(round(eigE_fpw[0], 2)), c="#B7094C")
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 1], label=r"$\psi_2$, $E_2=$" + str(round(eigE_fpw[1], 2)), c="#892B64")
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 2], label=r"$\psi_3$, $E_3=$" + str(round(eigE_fpw[2], 2)), c="#723C70")
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 3], label=r"$\psi_4$, $E_4=$" + str(round(eigE_fpw[3], 2)), c="#5C4D7D")
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 4], label=r"$\psi_5$, $E_5=$" + str(round(eigE_fpw[4], 2)), c="#2E6F95")
# plt.plot(x_fpw, 50 * eigPsi_fpw[:, 5], label=r"$\psi_6$, $E_6=$" + str(round(eigE_fpw[5], 2)), c="#0091AD")
# plt.plot(x_fpw, V_fpw, c="#151729")
#
# plt.title("Vezana stanja končne pot. jame globine 100 a.u z FD")
# plt.xlabel("x")
# plt.ylabel(r"$\psi_n(x)$")
# plt.legend(loc="upper right")
# plt.xlim(-2, 2)
# plt.show()

# Fuck it LHO
# lho_x, lho_num, lho_E, lho_V = lho(psi_init, upper, h_)
# n_solutions = np.shape(lho_E)[0]
# print("Found {} solutions!".format(n_solutions))
# # fpw_num = np.column_stack(fpw_num)
# plt.plot(lho_x, 5 * lho_num[0], label=r"$\psi_1$, $E_1=$" + str(round(lho_E[0], 2)), c="#B7094C")
# plt.plot(lho_x, 5 * lho_num[1], label=r"$\psi_2$, $E_2=$" + str(round(lho_E[1], 2)), c="#892B64")
# plt.plot(lho_x, 5 * lho_num[2], label=r"$\psi_3$, $E_3=$" + str(round(lho_E[2], 2)), c="#723C70")
# plt.plot(lho_x, 5 * lho_num[3], label=r"$\psi_4$, $E_4=$" + str(round(lho_E[3], 2)), c="#5C4D7D")
# plt.plot(lho_x, 5 * lho_num[4], label=r"$\psi_5$, $E_5=$" + str(round(lho_E[4], 2)), c="#2E6F95")
# plt.plot(lho_x, 5 * lho_num[5], label=r"$\psi_6$, $E_6=$" + str(round(lho_E[5], 2)), c="#0091AD")
# plt.plot(lho_x, lho_V, c="#151729")
# plt.title("Vezana stanja linearnega harmonskega oscilatorja")
# plt.xlabel("x")
# plt.ylabel(r"$\psi_n(x)$")
# plt.ylim(-10, 30)
# plt.xlim(-5, 5)
# plt.legend()
# plt.show()
# Fuck it LHO MK2
# lho_x, lho_num, lho_E, lho_V = fd_lho(1000, 3, 0)
# plt.plot(lho_num)
# plt.show()


# Dimension average error
fig, ax = plt.subplots()
stop = 1000
n_data, x_data, avg_err, img_sq = error_with_dim(2, stop, 0)
# plt.plot(n_data, avg_err)
# ani = ArtistAnimation(fig, img_sq, interval=30, repeat=True, blit=False)
# ani.save("psi1.mp4", "ffmpeg", fps=30)
# plt.title("Večanje natančnosti osnovnega stanja z večanjem")

# colors = cmr.take_cmap_colors('cmr.cosmic', stop, cmap_range=(0.15, 0.85), return_fmt='hex')
# plasma = plt.cm.get_cmap("plasma", stop)
#
#
# def update(iter):
#     plt.title("Osnovno stanje pri n = {}".format(iter))
#     if np.average(avg_err[iter]) > 0:  # Fix some orientations being down
#        return plt.plot(x_data[iter], avg_err[iter], c=plasma(iter))
#     else:
#        return plt.plot(x_data[iter], -avg_err[iter], c=plasma(iter))
#
# plt.xlabel("x")
# plt.ylabel(r"$\psi_1(x)$")
# ani = FuncAnimation(fig, update, frames=stop-3, interval=15, blit=False, repeat=False)
# plt.colorbar(plt.cm.ScalarMappable(cmap=plasma, norm=matplotlib.colors.Normalize(vmin=0, vmax=stop)), label="Dimenzija matrike")
# ani.save("psi0.mp4", "ffmpeg", fps=60)
# plt.show()


plt.plot(n_data, avg_err, c="#B7094C")
plt.title("Povprečno odstopanje od analitične rešitve")
plt.xlabel("n")
plt.ylabel("Avg. odstopanje")
plt.yscale("log")
plt.show()
