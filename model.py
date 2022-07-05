from types import SimpleNamespace
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import scipy.linalg as LA

# matrix functions:

def fill_diag(matrix):
    """Sets diagonal entries so each column sums to zero."""
    return matrix - np.diag(matrix.sum(axis=0))

def steady_state(matrix):
    """Compute steady state for a system y'=Ay"""
    w, v = LA.eig(matrix) # evecs & evals
    return np.abs(v[:, np.argmin(np.abs(w))]) # evec with eval closest to 0

# model setup

box_names = ["sed", "soil", "surf", "bio", "deep"]

# number of boxes
N = len(box_names)

# dictionary for readable indexing
box_dict = {}
for x in range(N):
    box_dict[box_names[x]] = x

#i.sed = 0, i.soil = 1, etc
i = SimpleNamespace(**box_dict)

# steady state (SS) magnitude vector (Tg P)
MAG_SS = np.zeros(N)
MAG_SS[i.sed] = 2e9
MAG_SS[i.soil] = 2e5
MAG_SS[i.surf] = 2800
MAG_SS[i.bio] = 44
MAG_SS[i.deep] = 1e5

GLOBAL_MAG_SS = MAG_SS.sum()

# initial magnitude vector (Tg P)
MAG_INIT = np.zeros(N)
MAG_INIT[i.sed] = 2e9
MAG_INIT[i.soil] = 2e5
MAG_INIT[i.surf] = 280
MAG_INIT[i.bio] = 44
MAG_INIT[i.deep] = 1e5

# scale initial conditions to match global magnitude
MAG_INIT *= (GLOBAL_MAG_SS / MAG_INIT.sum())

# time interval
T = np.arange(1e3) # yr

# biogeochemical parameters
WEATHERING_RATE = 2e4 # Tg sediment /yr
SED_ABUND = 0.001 # Tg P /Tg sediment
P_WEATHERING_RATE = WEATHERING_RATE * SED_ABUND # Tg P /yr
INSOLUBLE_FRAC = 0.9

SURF_GPP = 4e4 # Tg C /yr
C_MASS = 12 # g /mol
P_MASS = 31 # g /mol
REDFIELD_MOLE = 106
REDFIELD_MASS = REDFIELD_MOLE * C_MASS / P_MASS

REMIN_FRAC = .96

VERT_EXCH = 2 # m /yr
OCEAN_SA = 3.5e14 # m^2
VERT_FLOW = VERT_EXCH * OCEAN_SA / 1e12 # (1e12 g /Tg)

SURF_CONC = 0.025 # g /m^3
DEEP_CONC = 0.080 # g /m^3

# flux matrix
F = np.zeros((N, N))

# well constrained fluxes
# F[A,B] is from B to A
F[i.sed,i.soil] = P_WEATHERING_RATE * INSOLUBLE_FRAC
F[i.surf,i.soil] = P_WEATHERING_RATE * (1 - INSOLUBLE_FRAC)

F[i.bio,i.surf] = SURF_GPP / REDFIELD_MASS

F[i.surf,i.bio] = F[i.bio,i.surf] * REMIN_FRAC
F[i.deep,i.bio] = F[i.bio,i.surf] * (1 - REMIN_FRAC)

F[i.deep,i.surf] = VERT_FLOW * SURF_CONC
F[i.surf,i.deep] = VERT_FLOW * DEEP_CONC

# poorly constrained fluxes defined to satisfy SS for one box
F[i.sed,i.deep] = (F.sum(axis=1) - F.sum(axis=0))[i.deep]
F[i.soil,i.sed] = (F.sum(axis=1) - F.sum(axis=0))[i.sed]

# flux (Tg P /yr) to linear rate constants (/yr)
# divide each column of F (FROM) by source reservoir mag.
# and fill diagonal entries
K = fill_diag(F / MAG_SS)

# constant coeff. unforced system p' = Kp
def unforced_const_coeff(p, t):
    return K @ p

# solver
solution = odeint(unforced_const_coeff, MAG_INIT, T)

# ss solution for plot comparison
solution_ss = np.tile(MAG_SS, len(T)).reshape(len(T), N)

# plot
fig, ax = plt.subplots()
ax.plot(solution)
ax.set_yscale('log')
plt.legend(box_names, loc='upper right')
ax.plot(solution_ss, "--", linewidth=0.5, color='grey')
plt.savefig("figure.pdf")

# time-dependent forcing from Table 1 Avigad & Gvirtzman (2009)
erosion_forcing = {
    'time': np.array([635.0, 630.0, 615.0, 600.0, 530.0]), # Ma
    'crustal_thickness': np.array([50, 50, 42, 37, 35]), # km
    'mantle_lith_thickness': np.array([150, 0, 55, 85, 100]), # km
    'moho_temp': np.array([3, 1300, 900, 850, 650]), # deg C
    'crust_topo_contrib': np.array([-2.3, 3, 2.2, 1.7, 1.5]), # km
    'mantle_topo_contrib': np.array([0, 0, -0.5, -0.75, -1.3]), # km
    'calc_topo': np.array([0.5, 3, 1.7, 0.95, 0.2, ]), # km
}

erosion_forcing['time'] *= -1.0e9 # convert Ma to years

erosion_time_years = np.arange(erosion_forcing['time'][0], erosion_forcing['time'][-1], 100000)

interpolated_topo = interp1d(erosion_forcing['time'], erosion_forcing['calc_topo'])

def general(p, t):
    k_matrix = K
    b = np.zeros(N)
    return k_matrix @ p + b
