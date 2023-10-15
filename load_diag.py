import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as color
from src_fd import *

data = np.load("ver4_500.npz")
x = data["x"]
eigE = data["eigE"]
eigPsi = data["eigPsi"]

# plt.plot(x, get_err(eigPsi[:, 0], analytic(x, 1)), label=r"$\psi_1$", c="#B7094C")
# plt.plot(x, get_err(eigPsi[:, 1], analytic(x, 2)), label=r"$\psi_2$", c="#892B64")
# plt.plot(x, get_err(eigPsi[:, 2], analytic(x, 3)), label=r"$\psi_3$", c="#723C70")
# plt.plot(x, get_err(eigPsi[:, 3], analytic(x, 4)), label=r"$\psi_4$", c="#5C4D7D")
# plt.plot(x, get_err(eigPsi[:, 4], analytic(x, 5)), label=r"$\psi_5$", c="#2E6F95")
# plt.plot(x, get_err(eigPsi[:, 5], analytic(x, 6)), label=r"$\psi_6$", c="#0091AD")
# plt.title("Absolutne napake rešitev z FD in n = {}".format(x.shape[0]))
# plt.xlabel("x")
# plt.ylabel(r"$|num - analytic|$")
# plt.yscale("log")
# plt.legend()
# plt.show()

# Plot matrix
plt.imshow(eigPsi, cmap="plasma", norm=color.CenteredNorm(vcenter=0))
plt.title("Diagonalizirana matrika sistema n = {}".format(x.shape[0]))
plt.colorbar()
plt.show()
