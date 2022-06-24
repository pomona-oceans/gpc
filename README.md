# gpc

The purpose of this project is to explore the role of global phosphorous cycling in the unusual landscape of the late Neoproterozoic. We aim to account for a variety of factors including biogeochemistry, chemical erosion, and sedimentary deposition.

The basic construction is of the form:

$$\mathbf{y}' = A\mathbf{y} + \mathbf{b}.$$

In particular, $\mathbf{y}$ is a vector representing several reservoir magnitudes, so this is a system of equations. $A$ is matrix of constants and/or functions of time $t$, so the system is linear. $\mathbf{b}$ is a vector of forcing terms, which may be zero, (homogeneous) constants, or functions of $t$.

Each entry of $A$ and $\mathbf{b}$ will be set individually, since each flux and forcing term is specific to the nature of interconnections between phosphate pools like the surface ocean, sediments, biota, etc.
