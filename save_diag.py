import numpy as np
from src_fd import *

while True:
    dim = int(input("Matrix dimension: "))
    filename = str(input("Enter filename to save: "))
    check = str(input("Type \'yes\' to proceed: "))
    if check == "yes":
        break
    else:
        pass

x, eigE, eigPsi = fd_ipw(dim, 1, 0)
np.savez(filename, x=x, eigE=eigE, eigPsi=eigPsi)

