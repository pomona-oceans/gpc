import numpy as np

# reservoir names
boxes = ("sed", "soil", "surface_ocean", "ocean_biota", "deep_ocean")

N = len(boxes)

state_0 = np.zeros(N)

np.put(state_0, 0, 2000000000)

print(state_0)