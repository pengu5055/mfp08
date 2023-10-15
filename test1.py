import numpy as np
import matplotlib.pyplot as plt

n = 10
h = (5-0) / n

# Get A
A = np.zeros((n+1, n+1))
A[0, 0] = 1
A[n, n] = 1
for i in range(1, n):
    A[i, i-1] = 1
    A[i, i] = -2
    A[i, i+1] = 1

print(A)

# Get b
b = np.zeros(n+1)
b[1:-1] = -9.8*h**2
b[-1] = 50
print(b)

# solve the linear equations
y = np.linalg.solve(A, b)

t = np.linspace(0, 5, 11)

plt.figure(figsize=(10,8))
plt.plot(t, y)
plt.plot(5, 50, 'ro')
plt.xlabel('time (s)')
plt.ylabel('altitude (m)')
plt.show()

y_n1 = -9.8*h**2 + 2*y[0] - y[1]
(y[1] - y_n1) / (2*h)


def get_a_b(n):
    h = (np.pi / 2 - 0) / n
    x = np.linspace(0, np.pi / 2, n + 1)
    # Get A
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1
    A[n, n] = -2 + 4 * h ** 2
    A[n, n - 1] = 2
    for i in range(1, n):
        A[i, i - 1] = 1
        A[i, i] = -2 + 4 * h ** 2
        A[i, i + 1] = 1

    # Get b
    b = np.zeros(n + 1)
    for i in range(1, n + 1):
        b[i] = 4 * h ** 2 * x[i]

    return x, A, b


x = np.pi / 2
v = x - np.sin(2 * x)

n_s = []
errors = []

for n in range(3, 100, 5):
    x, A, b = get_a_b(n)
    y = np.linalg.solve(A, b)
    n_s.append(n)
    e = v - y[-1]
    errors.append(e)

plt.figure(figsize=(10, 8))
plt.plot(n_s, errors)
plt.yscale('log')
plt.xlabel('n gird points')
plt.ylabel('errors at x = $\pi/2$')
plt.show()