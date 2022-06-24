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

# reservoir names
boxes = ["sed", "soil", "surf", "biota", "deep"]
# number of boxes
N = len(boxes)

# for readable indexing
sed, soil, surf, biota, deep = np.arange(N)

# reservoir magnitute at steady state (Tg P)
reservoir_ss = np.array([2e9, 2e5, 2800, 44, 1e5])
total = reservoir_ss.sum()

# initialize flux matrix
P_flux = np.zeros((N,N))

# dependencies for flux definitions
weathering_rate = 2e4 # Tg P /yr
crustal_p_abundance = 0.001 # g P /g sediment
insol_frac = 0.9

redfield_C_P_molar = 106
P_molar_mass = 31 # g /mol
C_molar_mass = 12 # g /mol
redfield_C_P_mass = redfield_C_P_molar * C_molar_mass / P_molar_mass

gross_surf_prod = 4e4 # Tg C /yr

remin_frac = .96

vert_exchange = 2 # m /yr
surf_p_conc = 0.025 # g /m^3
deep_p_conc = 0.080 # g /m^3
ocean_surf_area = 3.5e14 # m^2
g_to_Tg = 1e12 # g /Tg

# define well constrained fluxes based on parameters
# entries set according to P_flux[TO,FROM]
P_flux[sed,soil] = weathering_rate * crustal_p_abundance * insol_frac
P_flux[surf,soil] = weathering_rate * crustal_p_abundance * (1 - insol_frac)
P_flux[biota,surf] = gross_surf_prod / redfield_C_P_mass
P_flux[surf,biota] = P_flux[biota,surf] * remin_frac
P_flux[deep,biota] = P_flux[biota,surf] * (1 - remin_frac)
P_flux[deep,surf] = vert_exchange * surf_p_conc * ocean_surf_area / g_to_Tg
P_flux[surf,deep] = vert_exchange * deep_p_conc * ocean_surf_area / g_to_Tg

# define poorly constrained fluxes based on steady state
P_flux[sed,deep] = (P_flux.sum(axis=1)-P_flux.sum(axis=0))[deep]
P_flux[soil,sed] = (P_flux.sum(axis=1)-P_flux.sum(axis=0))[sed]

# fill diagonal terms
fill_hollow(P_flux)

# fluxes at steady state (Tg P /y)
flux_ss = np.array([[ 0, 18,   0,   0,  2],
                    [20,  0,   0,   0,  0],
                    [ 0,  2,   0, 940, 56],
                    [ 0,  0, 980,   0,  0],
                    [ 0,  0,  18,  40,  0]])

# fill diagonal terms of flux_ss matrix
fill_hollow(flux_ss)

# flux (Tg P /yr) to linear rate constants (/yr)
K = P_flux / reservoir_ss

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
