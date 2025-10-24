import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from math import sqrt, pi
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

R = 15
r = 0.6


# original eigenvalues
def f(k, R, r):
    q1 = (k - pi * R**2) / (sqrt(2 * pi) * R)
    q2 = (k - pi * (r * R)**2) / (sqrt(2 * pi) * r * R)
    return 0.5 * (erfc(q1) - erfc(q2))


# supposed sorted version
def g(k, R, r):
    area = pi * R**2 * (1 - r**2)
    perimeter = 2 * pi * R * (1 + r)
    return 0.5 * erfc(sqrt(2 * pi) * (k - area) / perimeter)


# Generate a list of k values
k_values = np.array(range(int(1.4 * pi * R**2)))

# Calculate f(k), g(k) for each k
f_values = f(k_values, R, r)
g_values = g(k_values, R, r)

# Sort f(k) in decreasing order
sorted_f_values = np.sort(f_values)[::-1]


def k_for_level(level):
    """Return the k-value where g(k) crosses the supplied level."""
    if level > g_values.max() or level < g_values.min():
        return None
    # g_values is decreasing in k, so reverse to satisfy interp requirement
    return float(np.interp(level, g_values[::-1], k_values[::-1]))


k_low = k_for_level(0.55)
k_center = k_for_level(0.5)
k_high = k_for_level(0.45)

# Fallback to a fixed window if interpolation failed
if None in (k_low, k_center, k_high):
    k_low, k_center, k_high = 550, 600, 650

fig, ax = plt.subplots(figsize=(11, 3.0), dpi=150)

ax.plot(k_values, f_values, color='tab:blue', linewidth=1.5, label=r'Unsorted $f(k)$')
ax.plot(k_values, sorted_f_values, color='tab:purple', linewidth=2.5, label=r'Sorted $f(k)$ (descending)')
ax.plot(
    k_values,
    g_values,
    color='tab:orange',
    linewidth=2,
    linestyle='--',
    label=r'$0.5 \cdot \mathrm{erfc}\left(\frac{k - \pi R^2(1 - r^2)}{\sqrt{2\pi} R(1+r)}\right)$',
)

ax.axvspan(k_low, k_high, color='0.85', alpha=0.5, lw=0)
ax.set_title(f'Comparison of $f(k)$ Variants (r = {r}, R = {R})')
ax.set_xlabel('$k$')
ax.set_ylabel('$f(k)$')
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.9))

# Zoomed inset centered around where g(k) hits 0.5
zoom_ax = inset_axes(ax, width='42%', height='42%', loc='upper right', borderpad=1.0)
# Nudge the inset a little to the right without raising it further from the top edge.
inset_pos = zoom_ax.get_position()
zoom_ax.set_position([inset_pos.x0 + 0.015, inset_pos.y0, inset_pos.width, inset_pos.height])
zoom_ax.set_box_aspect(1)
zoom_ax.plot(k_values, f_values, color='tab:blue', linewidth=1.5)
zoom_ax.plot(k_values, sorted_f_values, color='tab:purple', linewidth=2.5)
zoom_ax.plot(
    k_values,
    g_values,
    color='tab:orange',
    linewidth=2,
    linestyle='--',
)

zoom_ax.set_xlim(k_low, k_high)
zoom_ax.set_ylim(0.45, 0.55)
zoom_ax.grid(True, linestyle=':', alpha=0.6)
zoom_ax.set_xticks([round(k_low), round(k_center), round(k_high)])
zoom_ax.set_yticks([0.45, 0.5, 0.55])
zoom_ax.tick_params(labelsize=8)

mark_inset(ax, zoom_ax, loc1=2, loc2=4, fc='none', ec='0.4')

fig.tight_layout()

fig.savefig('sorted.png', bbox_inches='tight')
#plt.show()
