import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("shoot_try.csv")
print(data)
plt.plot(data)
plt.show()
