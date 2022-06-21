import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# functions

def row_col_test(matrix):
    """True if all columns and rows sum to zero."""
    return matrix.sum(axis=0) == matrix.sum(axis=1) == 0

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
flux_ss -= np.diag(flux_ss.sum(axis=0))

# convert fluxes to linear rate constants
K = flux_ss / reservoir_ss

def f(y, t):
    return K @ y

t_int = np.arange(5)

solution = odeint(f, reservoir_ss, t_int)

fig, ax = plt.subplots()

ax.plot(t_int, solution,)

ax.set_yscale('log')

plt.savefig("figure.pdf")
