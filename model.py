import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# functions:

def fill_hollow(matrix):
    """Sets diagonal entries so each column sums to zero."""
    matrix -= np.diag(matrix.sum(axis=0))

def row_col_test(matrix):
    """True if all columns and rows sum to zero."""
    return matrix.sum(axis=0) == matrix.sum(axis=1) == 0

# - - - - - - - - - - - - - -

# characterize model:

# reservoir names
boxes = ["sed", "soil", "surface_ocean", "ocean_biota", "deep_ocean"]
N = len(boxes)

# reservoir magnitute at steady state (Tg P)
reservoir_ss = np.array([2e9, 2e5, 2800, 44, 1e5])
total = reservoir_ss.sum()

# fluxes at steady state (Tg P /y)
flux_ss = np.array([[ 0, 18,   0,   0,  2],
                    [20,  0,   0,   0,  0],
                    [ 0,  2,   0, 940, 56],
                    [ 0,  0, 980,   0,  0],
                    [ 0,  0,  18,  40,  0]])

# fill diagonal terms of flux_ss matrix
fill_hollow(flux_ss)

# flux (Tg P /yr) to linear rate constants (/yr)
K = flux_ss / reservoir_ss

# constant coeff. linear homogeneous system y' = Ky
def f(y, t):
    return K @ y

# initial reservoir conditions
reservoir_init = np.array([2e9, 2e5, 280, 44, 1e5]) # Tg P

# scale to match global magnitude
reservoir_init *= (total / reservoir_init.sum())

# time interval
t_int = np.arange(1e3) # yr

# solver
solution = odeint(f, reservoir_init, t_int)

# ss solution for plot comparison
solution_ss = np.tile(reservoir_ss, len(t_int)).reshape(len(t_int), N)

# plot
fig, ax = plt.subplots()
ax.plot(solution)
ax.set_yscale('log')
plt.legend(boxes, loc='upper right')
ax.plot(solution_ss, "--", linewidth=0.5, color='grey')
plt.savefig("figure.pdf")
