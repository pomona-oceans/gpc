from types import SimpleNamespace as SN
import numpy as np
import pandas as pd
from scipy.integrate import odeint
# from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# boxes
boxes = pd.Series(["sed", "soil", "surf", "bio", "deep"])
N = len(boxes)

# readable indexing: i.sed returns 0, etc.
i = SN(**dict((v, k) for k, v in boxes.items()))

# Tg P at steady state (SS)
P_SS = np.zeros(N)
P_SS[i.sed] = 2e9
P_SS[i.soil] = 2e5
P_SS[i.surf] = 2800
P_SS[i.bio] = 44
P_SS[i.deep] = 1e5

# Global Tg P
GLOBAL_P = P_SS.sum()

# Tg P when t = 0 (initial conditions)
P_0 = np.zeros(N)
P_0[i.sed] = 2e9
P_0[i.soil] = 2e6
P_0[i.surf] = 2800
P_0[i.bio] = 44
P_0[i.deep] = 1e5

# scale initial conditions to match global reservoir
P_0 *= (GLOBAL_P / P_0.sum())

# time interval
t = np.arange(start = 0, stop = 1e4, step = 1) # yr

# biogeochemical parameters
WEATHERING_RATE = 2e4 # Tg sediment /yr
SED_ABUND = 0.001 # Tg P /Tg sediment
P_WEATHERING_RATE = WEATHERING_RATE * SED_ABUND # Tg P /yr
INSOLUBLE_FRAC = 0.9

SURF_GPP = 4e4 # Tg C /yr
C_MASS = 12 # g /mol
P_MASS = 31 # g /mol
REDFIELD_MOLE = 106 # mol /mol
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
# F[i.A, i.B] is from B to A
F[i.sed, i.soil] = P_WEATHERING_RATE * INSOLUBLE_FRAC
F[i.surf, i.soil] = P_WEATHERING_RATE * (1 - INSOLUBLE_FRAC)

F[i.bio, i.surf] = SURF_GPP / REDFIELD_MASS

F[i.surf, i.bio] = F[i.bio, i.surf] * REMIN_FRAC
F[i.deep, i.bio] = F[i.bio, i.surf] * (1 - REMIN_FRAC)

F[i.deep, i.surf] = VERT_FLOW * SURF_CONC
F[i.surf, i.deep] = VERT_FLOW * DEEP_CONC

# poorly constrained fluxes defined to satisfy one SS condition
# axis=1: influx (row); axis=0: outflux (col.)
F[i.sed, i.deep] = (F.sum(axis = 1) - F.sum(axis = 0))[i.deep]
F[i.soil, i.sed] = (F.sum(axis = 1) - F.sum(axis = 0))[i.sed]

# rate constant matrix (/yr)
K = F / P_SS

# set diagonal terms so each col. sums to 0
K -= np.diag(K.sum(axis=0))

# constant coeff. unforced system p' = Kp
def unforced_const_coeff(p, t):
    return K @ p

# solver
solution = odeint(unforced_const_coeff, P_0, t)

# plot
fig, ax = plt.subplots()
ax.plot(t, solution)
ax.plot(t[-1], [P_SS], 'x')

ax.set_yscale('log')
# ax.set_xscale('log')

plt.legend(boxes)

plt.savefig("figure.pdf")
