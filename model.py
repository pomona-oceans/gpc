import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# reservoir names
boxes = ("sed", "soil", "surface_ocean", "ocean_biota", "deep_ocean")
N = len(boxes)

# reservoir magnitute at steady state (Tg P)
reservoir_ss = np.array([2000000000, 200000, 2800, 44, 100000])

# fluxes at steady state (Tg P / y)
flux_ss = np.array([[ 0, 18,   0,   0,  2],
                    [20,  0,   0,   0,  0],
                    [ 0,  2,   0, 940, 56],
                    [ 0,  0, 980,   0,  0],
                    [ 0,  0,  18,  40,  0]])

# fill diagonal terms of flux matrix
flux_ss += np.diag(flux_ss.sum(axis=0))

# convert fluxes to linear rate constants
K = (flux_ss.T / reservoir_ss).T

def ss_test(matrix):
    """True if all columns and rows sum to zero."""
    return matrix.sum(axis=0) == matrix.sum(axis=1) == 0

def f(y, t):
    return K @ y

t = np.arange(100)

solution = odeint(f, reservoir_ss, t)

fig, ax = plt.subplots()

ax.plot(t, solution, linewidth=2.0)

plt.savefig("figure.pdf")
