This folder contains the code necessary to run the solver and process the output:

- nondimensionalse (called by parameters.m) converts the user-input parameters to dimensionless quantities
- construct_protocol.m (called by parameters.m) converts a set of instructions for the illumination and voltage protocols into the dimensionless functions of time `light` and `psi` required by the solver and creates the dimensionless `time` and `splits` vectors
- numericalsolver.m (called by master.m) calls FE_solve.m, calculate_currents.m and then, if it exists, completion_tasks.m before re-dimensionalising and packaging up the solution
- FE_solve.m (called by numericalsolver.m) performs the solution procedure, using the scripts contained in Code/Solver
- calculate_currents.m (called by numericalsolver.m) is used to calculate the total dimensionless current density as well as the losses due to interfacial recombination
- find_Voc.m and precondition.m are functions (called by FE_solve.m if needed) to obtain alternative initial conditions when the voltage protocol starts from either open-circuit or a fixed voltage not equal to the cell's built-in voltage, respectively
- create_stats_funcs.m (called by nondimensionalise.m) is a function to create a statistical model and provide functions to quickly and accurately evaluate statistical functions and their inverses
- FourierFit.m (called by impedance_analysis.m) is a function that fits a sinusoidal curve to some periodic input signal using the Fourier transform
- impedance_analysis.m (called by IS_solver.m) is a function that calculates the impedance from a structure array of impedance measurement simulations
- IS_solver.m (called by master.m) is the intermediary step between master.m and numericalsolver.m for impedance simulations, iteratively replacing the voltage protocol with a sinusoidal protocol and calling numericalsolver.m to obtain the solution for each sample frequency

In addition:
- struct2array.m is used (in almost all functions) to extract the contents of a structure (such as `params` or `sol`)
- PrintVolt.m and PrintTime.m are output functions for ode15s, which allow ode15s to output its progress throughout the solver procedure
- sec2hms.m is used to convert the run time in seconds to hours, minutes and seconds
- completion_tasks.m can be used to run additional commands after the solution is found
