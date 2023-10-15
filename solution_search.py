import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("shoot_try2.csv", delimiter=",")
x = np.genfromtxt("E_arr.csv")
eigE = [9.8696044, 39.4784176, 88.82643961, 157.91367042, 246.74011003, 355.30575844, 483.61061565]
eigE_round = [round(element, 2) for element in eigE]
zero = [0, 0, 0, 0, 0, 0, 0]
fig, ax = plt.subplots()
plt.plot(x, data, c="#2E6F95", label="Ciljanje robnega pogoja")
plt.scatter(eigE, zero, c="#B7094C", label="Rešitve")
ax.annotate(eigE_round[0], (eigE[0], 0), xytext=(eigE[0] - 5, .3), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[1], (eigE[1], 0), xytext=(eigE[1] - 5, .2), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[2], (eigE[2], 0), xytext=(eigE[2] - 10, .3), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[3], (eigE[3], 0), xytext=(eigE[3] - 10, .2), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[4], (eigE[4], 0), xytext=(eigE[4] - 10, .2), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[5], (eigE[5], 0), xytext=(eigE[5] - 30, .2), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))
ax.annotate(eigE_round[6], (eigE[6], 0), xytext=(eigE[6] - 30, .2), bbox=dict(boxstyle="round", fc="#6ed6eb"), arrowprops=dict(arrowstyle="-", Connectionstyle="angle,angleA=0,angleB=-90,rad=10"))

plt.title("Iskanje rešitev z strelsko metodo")
plt.xlabel("E [a.u.]")
plt.ylabel("Odstopanje od robnega pogoja")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.legend()
plt.show()

# plt.plot(data[:, 0] + data[:, 1])
# plt.plot(data[:, 0] - data[:, 1])
# plt.axhline(alpha=1, ls=":", c="#adadad")
# plt.show()
