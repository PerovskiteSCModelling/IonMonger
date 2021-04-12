# IonMonger

A drift-diffusion model for ion migration and charge carrier transport across a planar perovskite solar cell (PSC).

This code can be used to simulate the internal state of a PSC over time. The three core layers of a PSC, namely the electron transport layer, perovskite absorber layer and hole transport layer, are modelled explicitly in one spatial dimension. The model variables are the electric potential, halide ion vacancies (existing only within the perovskite layer), electrons (within the ETL and perovskite layers) and holes (within the perovskite and HTL). A variety of experimental protocols can be simulated, including changes in the applied voltage and/or illumination intensity that occur over timescales on the order of microseconds to minutes. The code also outputs the current density and voltage which can be used to plot the current-voltage characteristics of a PSC, including current-voltage hysteresis due to the movement of halide ion vacancies. Please read the [GUIDE](GUIDE.md) to get started.


# Use Cases

This code is intended for use by researchers in the field of perovskite solar cells. Example use cases include:

- simulating current-voltage curves, with the ability to change key material properties in order to investigate trends in performance and the extent of hysteresis
- simulating photo-current or photo-voltage transients to investigate the effects of halide ion migration
- visualising the effects of halide ion migration on the internal electrical state of a PSC

The authors of this code published an investigation into how material properties of the transport layers affect perovskite solar cell performance in [Energy & Environmental Science](https://doi.org/10.1039/C8EE01576G), while working at the Universities of Southampton, Bath and Portsmouth.


# Requirements and Other Information

Requirements: MATLAB (version R2018b).

This code was first created at the University of Southampton in 2016. See [AUTHORS](AUTHORS.md) for a list of contributors and [LICENSE](LICENSE) for the conditions of use.

For details of the changes to the code since the first release, see the Changelog on the (IonMonger Wiki)[https://github.com/PerovskiteSCModelling/IonMonger/wiki].

If you encounter a problem or any unexpected results, please create an Issue on the GitHub website, add details of the problem (including the error message and MATLAB version number) and attach the parameters.m file in use when the problem occurred. For other enquiries, please contact N.E.Courtier(at)soton.ac.uk.


# How to Cite this Code

Please cite the release paper published in the [Journal of Computational Electronics](https://link.springer.com/article/10.1007/s10825-019-01396-2) by using the [citation.bib](citation.bib) file.


# Technical Features

This code is based on the finite element scheme first presented in our paper in [Applied Mathematical Modelling](https://doi.org/10.1016/j.apm.2018.06.051) and is performed on a non-uniform ("tanh") spatial grid.

Files in the main folder:
  - master.m for running a single simulation
  - parameters.m for setting the inputs to the simulation
  - reset_path.m adds all subfunctions to the MATLAB path

The Code/ folder contains all subfunctions, including
  - a function to turn a list of instructions into a protocol for the light or applied voltage
  - functions that provide the ability to find the steady-state Voc and simulate open-circuit conditions
  - a function to plot current-voltage (`J`-`V`) data as well as the recombination currents (`Jl`, `Jr`)

The solution is saved in dimensional form into one output file, which also contains the input data.
