import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from math import sqrt, pi

R = 20
r = 0.7

#original eigenvalues
def f(k, R, r):
    q1 = (k - pi * R**2) / (sqrt(2 * pi) * R)
    q2 = (k - pi * (r * R)**2) / (sqrt(2 * pi) * r * R)
    return 0.5 * (erfc(q1) - erfc(q2))

#supposed sorted version
def g(k, R, r):
    area = pi * R**2 * (1 - r**2)
    perimeter = 2 * pi * R * (1 + r)
    return 0.5 * erfc(sqrt(2*pi)*(k -area)/perimeter)

# Generate a list of k values
k_values = np.array(range(int(1.4 * pi * R**2)))

# Calculate f(k), g(k) for each k
f_values = f(k_values, R, r)
g_values = g(k_values, R, r)

# Sort f(k) in decreasing order
sorted_f_values = np.sort(f_values)[::-1]

plt.figure(figsize=(12, 5), dpi=300)

# Plot 1: Original f(k) vs k
plt.subplot(2, 1, 1)
plt.plot(k_values, f_values, color='blue', label='$f(k)$ vs $k$')
plt.title(f'Unsorted $f(k)$, r = {r}, R = {R}')
plt.xlabel('$k$')
plt.legend()
plt.grid(True)

# Plot 2: Sorted f(k) and Standalone erfc on the same axes
plt.subplot(2, 1, 2)
plt.plot(sorted_f_values, color='green', label='Sorted $f(k)$ (Descending)')
plt.plot(g_values, color='red', linestyle='--', label=r'Standalone $0.5 \cdot \operatorname{erfc}\left(\frac{k - \pi R^2(1 - r^2)}{\sqrt{2\pi} R(1+r)}\right)$')
plt.title('Sorted $f(k)$ and $\operatorname{erfc}$ from proposition')
plt.xlabel('$k$')
plt.legend()
plt.grid(True)

plt.tight_layout()

plt.show()
