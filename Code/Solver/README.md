
This folder contains the following files:
- create_matrices.m (takes `params`, `vectors`; returns `matrices`)
- create_vectors.m (takes `params`; returns `vectors`)
- Jac.m (takes `params`; returns Jacobian matrix)
- mass_matrix.m (takes `params`; returns mass matrix)
- RHS.m (takes `time`, state vector, potential profile, `params`, `vectors`, `matrices`; returns derivative of state vector wrt time)
- initial_conditions.m (takes `vectors`, `params`, `matrices`; returns the initial conditions vector `sol_init`)
- apply_Poisson.m (takes `sol_init`, `vectors`, `params`, `matrices`; returns `sol_init`)
- unpack.m (takes `numsol`, `params`; returns the solution structure `dstrbns`)

Jac.m, mass_matrix.m and RHS.m can also accept an additional optional input argument which, when set equal to `'open-circuit'`, modifies them to solve the system of equations at open-circuit (Voc), rather than a fixed applied voltage.

Additionally an optional completion_tasks.m file (which can be found in Code/Common) may be specified, which will be run at the end of the solve. Use this to output ion mass conservation messages etc.

The finite element method is detailed in our paper in [Applied Mathematical Modelling](https://doi.org/10.1016/j.apm.2018.06.051) and in the release paper published in the [Journal of Computational Electronics](https://link.springer.com/journal/10825).