# User Guide for IonMonger

This guide explains how to use IonMonger by editing a parameters.m file and running simulations and analysis from the command window. Users without any MATLAB knowledge can alternatively use IonMonger Lite, a live scipt that serves as an interface to run current-voltage sweep and impedance spectroscopy simulations without the need to edit raw code. Note that IonMonger Lite cannot offer the flexibility or functionality of the full version.

## Contents

- Quick Start Guide

- Example Solutions

- Details of Model

- How to Edit the Input File

- How to Analyse the Ouput

- Common Errors

- How to Cite this Code


# Quick Start Guide

This code can be used to simulate the internal state, current-voltage behaviour, and impedance spectra of a three-layer planar perovskite solar cell. The computed variables are the ion vacancy density, the electric potential and the charge carrier (electron and hole) concentrations.

An example solution may be obtained and plotted using the following procedure.

1) Download the git repository and make the folder containing the [master.m](master.m) file the current folder in Matlab.

2) Make a copy of the [parameters_template.m](parameters_template.m) file and rename it parameters.m.

3) Enter `master;` in the command window. The program should now run using the default set of inputs given in parameters.m, which includes the model parameters, simulation settings and experimental protocol. Notes on amending these inputs are given below.

The command window will show the progress of the solver during the calculation and when it is finished. This example simulation should take less than a minute to complete. Automatically-generated figures 98 and 99 show the light and voltage protocols that are being simulated. To suppress these outputs, change the `Verbose` setting in parameters.m to `false`.

4) To plot the solution, one can use or adapt any of the following options:

a) To plot the example output (a J-V curve), enter `plot_sections(sol,[2,3]);` in the command window.

The new figure shows the photo-generated current density (`J`) versus applied voltage (`V`) in orange, with a solid line showing the reverse scan and a dashed line showing the forward scan. Additionally, the losses due to interfacial recombination at the ETL/perovskite interface (`Rl` in blue) and the perovskite/HTL interface (`Rr` in red) are shown, with dot-dashed lines for the reverse scans and dotted lines for the forward scans.

b) To plot the transient current density, enter the following statement in the command window.

```
figure; plot(sol.time,sol.J);
xlabel('Time (s)'); ylabel('Current density (mA/cm2)');
```

c) To plot the evolution of any of the solution variables across the perovskite layer (the ion vacancy density `P`, electric potential `phi`, electron concentration `n`, hole concentration `p`), enter e.g.

```
figure; surf(sol.vectors.x, sol.time, sol.dstrbns.phi,'LineStyle','none');
xlabel('Distance (nm)'); ylabel('Time (s)'); zlabel('Electric potential (V)');
```

in the command window. Similarly, to plot either of the solution variables across the ETL (electric potential `phiE`, electron concentration `nE`) or across the HTL (electric potential `phiH`, hole concentration `pH`), replace `x` with either `xE` or `xH`, respectively.

d) To generate more advanced plots, the whole solution can be loaded into the Matlab environment using:

```
reset_path; load('Data/simulation.mat');
```

The solution structure `sol` contains both the input (`params`, `vectors`) and output (`dstrbns`, `J`). To see the contents of a structure, enter e.g. `sol` or `sol.dstrbns` in the command window.

e) To animate the solution and create an MP4 file, use the function `animate_sections`. This requires the arguments `sol`, a solution structure; `sections`, a vector of integers specifying which sections of the protocol to animate; `length`, the video length in seconds and `filename`, a string under which the video will be saved. Additional optional arguments can be given that determine the frame rate and resolution of the video. The defaults are 30fps and 1920x1080p. For example:
```
animate_sections(sol,[2,3],10,'Videos\example_video')
```

5) To run further simulations, open the parameters.m file and edit as described in the following section entitled How to Edit the Input File. Then repeat the procedure from step 3.


# Example Solutions

Three example solutions are included in the release paper published in the [Journal of Computational Electronics](https://link.springer.com/article/10.1007/s10825-019-01396-2). To view more results produced using an equivalent code, please refer to our paper in [Energy & Environmental Science](https://doi.org/10.1039/C8EE01576G).


# Details of Model

The model is described in the release paper in the [Journal of Computational Electronics](https://link.springer.com/article/10.1007/s10825-019-01396-2). A detailed analysis of the performance of different solution methods, in terms of speed and numerical accuracy, is given for a comparable, single-layer model in our paper in [Applied Mathematical Modelling](https://doi.org/10.1016/j.apm.2018.06.051).


# How to Edit the Input File

## How to Change the Simulation Settings/Parameters

All simulation settings and model parameters are listed in the [parameters.m file](parameters_template.m file). To change a setting or a value of a parameter, edit and save the file as parameters.m. The code is not guaranteed to obtain a solution for all combinations of parameter values and simulations protocols, please refer to the section on Common Errors for advice on avoiding known causes of failure.


## How to Construct a Simulation Protocol

Here, the term 'simulation protocol' refers to the experimental conditions being applied to the cell. The solver requires two functions of time, namely `light` and `psi`, to describe the illumination intensity and the applied voltage during the measurement. The two functions must be defined for the duration of the simulation, described by the vector of time points `time`.

A quick way to define any simple protocol is to use the construct_protocol.m function. This function converts a set of instructions into continuous functions for `light` and `psi` and also outputs a suitable `time` vector.
The function can be used to generate the necessary solver input to describe any protocol that can be split into a finite number of sections, each of which has either a linear or tanh/cosine-like form. The tanh option can be used to describe a smooth increase/decrease that starts off quickly before flattening off (in the shape of tanh(t), t=0..3), while the cosine option creates a sigmoidal transition between two values (in the shape of cos(t), t=-pi..pi).

Either a single value for the initial applied voltage or `'open-circuit'` must be included as the first entry in the list of instructions for the applied voltage. The set of instructions required for each section of the protocol that follows are: the type of section (`'linear'`/`'tanh'`/`'cosine'`), the duration of the section in seconds and the value of the applied voltage at the end of the section, i.e.

{'type of slope', duration in seconds, voltage to end at}

For example, the set of instructions to describe a voltage protocol that starts from the cell's built-in voltage and consists of an initial 5-second preconditioning step in which the applied voltage is smoothly increased to 1.2 V, followed by a 100 mV/s J-V measurement (a reverse scan to short-circuit, immediately followed by a forward scan) is:

`{Vbi, 'tanh', 5, 1.2, 'linear', 1.2/0.1, 0, 'linear', 1.2/0.1, 1.2}`

Note that the scan rate can be converted to a duration using the following relation:
Duration (s) = Voltage Change (V)/Scan Rate (V/s).

A function to describe a varying illumination protocol can be constructed in a similar way, including a single value for the initial illumination intensity as the first entry in the list of instructions. A value of 1 corresponds to the light intensity that provides the flux of photons (with energy above the band gap) entering the perovskite layer defined by `Fph` in the parameters file. As an example, the set of instructions for a 0.1-second pulse of light after 10 seconds in the dark followed by another 10 seconds in the dark is:

`{0, 'linear', 10, 0, 'linear', 0.1, 1, 'linear', 10, 0}`

For a constant `light_intensity` or `applied_voltage`, only the starting value needs to be entered (however at least one of the protocols must describe a length of time).

The automatically generated `time` vector is composed of 100 time steps per section such that the first section is described by time(1:101), the second is described by time(101:201) and so on. Within each section, the time points can be either uniformly or logarithmically spaced, by setting `time_spacing` equal to either `'lin'` or `'log'`, respectively, in parameters.m.

In addition to the `time` vector, another vector of time points, namely `splits`, is automatically generated by construct_protocol.m. This vector includes all the times which mark the end or beginning of a section. In order to cope with possible abrupt changes in the gradient of the input functions, the solver makes separate calls to ode15s for each section. The [Matlab documentation](https://mathworks.com/help/matlab/ref/odeset.htmlmathworks.com/help/matlab/ref/ode15s.html) for ode15s advises the user to take such an approach. To avoid making separate calls to ode15s and instead attempt to calculate the solution in one call, set the `UseSplits` option in parameters.m to `false`.

The last output, `findVoc`, is automatically set to `true` if the applied voltage begins with `'open-circuit'` and `false` otherwise. When this option is enabled, the solver will attempt to locate the steady-state value of Voc, then adjust the first section of the voltage protocol automatically before continuing with the standard solution procedure.

To describe a more complex experimental protocol, the user can generate their own `light` and `psi` functions and corresponding vectors of time points `time` and `splits`. Note that these four inputs must be dimensionless and must be included as fields in the `params` structure in order to be passed to the solver.



## How to Model Open-Circuit Conditions

It is also possible to use this code to simulate open-circuit conditions. To simulate a cell that is kept at open-circuit throughout a measurement, the value of the `applied_voltage` (if using construct_protocol.m) must be set equal to `{'open-circuit'}`. In this case, a protocol for the light intensity must be specified, e.g. one can use `{1, 'linear', 1, 1}` for constant illumination for 1 second.

The only mathematical modification required to model open-circuit conditions is to replace the boundary conditions on the electric potential at the metal contacts (all other equations and boundary conditions remain the same) as described in the release paper. These modifications are contained within four scripts in the Code/Solver folder, namely Jac.m, AnJac.m, mass_matrix.m and RHS.m.


## How to Model Impedance Spectroscopy Measurements

Users can perform impedance spectroscopy simulations by changing the voltage protocol in the parameters file. If the first entry in `applied_voltage` is set to `'impedance'`, `master.m` will call `IS_solver.m`, which simulates the cell in steady state at the DC voltage, and then automatically creates the sinusoidal voltage protocol for each sample frequency and calls `numericalsolver.m` to obtain the solution. `applied_voltage` requires several more entries (all floating point numbers) to specify the bounds of the impedance measurements. These are, in order, the smallest sample frequency (Hz); the largest sample frequency (Hz); the DC applied voltage (V); the AC applied voltage amplitude (V); the number of frequencies to sample (must be a non-negative integer); the number of complete sine waves to apply at each frequency.
For example, to make 70 impedance measurements at frequencies between 1mHz and 1MHz, with a DC voltage of 0.8V, AC perturbations of amplitude 20mV and performing 6 complete sine waves at each frequency:
```
applied_voltage = {'impedance', ...
    1e-4, ...   % minimum impedance frequency (Hz)
    1e7, ...    % maximum impedance frequency (Hz)
    0.6, ...    % DC voltage (V)
    20e-3, ...  % AC voltage amplitude (V)
    6, ...      % number of frequencies to sample
    5, ...      % number of sine waves
    10};        % time to hold at DC voltage (s)
```
The sample frequencies will be logarithmically spaced between the minimum and maximum frequencies.

It is recommended that users maintain a constant light intensity throughout impedance measurements. For example, `light_intensity = {1}` for 1 Sun equivalent light intensity throughout.

For more information about how to analyse impedance spectra, see the section titled 'How to Analyse the Ouput'.


## How to Change the Statistical Models of Carriers in the Transport Layers

Users can choose to employ degenerate statistical models to carriers in the transport layers, where carrier densities often become high enough that the Boltzmann approximation is innaccurate. Low carrier concentrations in the perovskite mean that the Boltzmann approximation is considered to be accurate here. The `stats` structure in the parameters file contains all the information necessary to construct the statistical models. This structure has two fields, `ETL` and `HTL`. If either of these fields are missing, the statistical model in that region will be automatically set to be the Boltzmann approximation to Fermi-Dirac statistics and the transport model will be unchanged from v1.4 of IonMonger. To apply a specific statistical model in one of the transport layers, the corresponding field of `stats` must have the field `model` and the field `Boltzmann`. `model` determines which statistical model should be employed. The currently supported models are `'FermiDirac'`, `'GaussFermi'`, and `'Blakemore'`. `Boltzmann` can be set to `true` to use the Boltzmann approximation to the statistical model in question or `false` to use the full degenerate model. The Gauss-Fermi model also requires the field `s`, giving the level of Gaussian disorder in units of eV. The Blakemore model requires the additional field `lim` giving the (dimensionless) limiting concentration. If the Blakemore statistical function is written as B(x) = 1/(exp(-x) + gamma), then lim=1/gamma.


# How to Analyse the Ouput

The solution is saved in dimensional form into one output file. This output .mat file takes the form of a structure called `sol`. There are three possible forms of the `sol` structure, depending on the simulation protocol.

1) For transient simulations, `sol` contains:
  - the `vectors` structure, containing column vectors `x`, `dx`, `xE`, `dxE`, `xH` and `dxH`
  - the `params` structure, containing all input and calculated parameters
  - the `dstrbns` structure, containing `P`, `phi`, `n`, `p`, `phiE`, `nE`, `phiH` and `pH` (note that each variable is stored in the form P(t,x) and the command `surf(sol.vectors.x,sol.time,sol.dstrbns.P);` plots the whole evolution)
  - the `time`, `V`, `J`, `Jl`, `Jr`, `Jd` vectors of the same length
  - the `timetaken` for the simulation

2) For impedance spectroscopy simulations, `sol` is a structure array, where each structure is the full solution structure for a single frequency measurement, with the additional field `impedance_protocol`, the cell containing the parameters of the impedance protocol as `params.applied_voltage` will be replaced for each frequency with a sinusoidal voltage protocol. The frequency of the n-th simulation can be found by `1/sol(n).params.applied_voltage{2}`.

When taking a large number of sample frequencies, the solution structure array can become a large file, with much of the data being unnecessary.

3) To reduce the size of the solution, set `reduced_output` to `true` in the parameters file. In this case, `IS_solver.m` will perform the impedance analysis and discard most of the data before returning the solution. In this case, `sol` is a single structure that contains:
  - the `vectors` structure, containing column vectors `x`, `dx`, `xE`, `dxE`, `xH` and `dxH`
  - the `params` structure, containing all input and calculated parameters
  - the `dstrbns` structure, containing the steady state solutions of `P`, `phi`, `n`, `p`, `phiE`, `nE`, `phiH` and `pH` at the DC voltage. (note that each variable is stored in the form P(x) and the command `plot(sol.vectors.x,sol.dstrbns.P);` plots the distribution)
  - the `J`, `Jl`, `Jr` values in steady state at the DC voltage.
  - the vectors `freqs`, `R`, and `X`, containing the frequency, and the real and imaginary parts of the impedance, respectively, for each sample frequency.

The user must then choose how to analyse or plot the data. One example plotting function, called plot_sections.m, can be found in the Code/Plotting folder. Any solution of an impedance spectroscopy can be given to `plot_IS.m` to generate Nyquist and Bode plots as well as any user-defined plots. In order to perform an identical analysis at the end of every simulation, the user can add commands to the completion_tasks.m function, which can be found in the Code/Common folder. FourierFit.m is a function that can be used to fit a sinusoidal curve to a periodic signal. Example code is included in this function to analyse the phase of any variable from an impedance simulation.


# Common Errors

We find that the most common errors returned by Matlab's built-in solver `ode15s` are that it has failed either at the start or at a particular time during the integration. Such errors return one of the following messages in orange or red text. Either

> Error using daeic3 (line 230)
Need a better guess y0 for consistent initial
conditions.

or

> Warning: Failure at t=.... Unable to meet integration tolerances without reducing the step size below the smallest value allowed (...) at time t.

> Output argument "dstrbns" (and maybe others) not assigned during call to "FE_solve".

These errors can sometimes be circumvented by adjusting the values of `N`, `atol` and `rtol`, which can be found in parameters.m. The parameter `N` controls the number of grid points on which the solution is computed, while `atol` and `rtol` are used to set the absolute and relative error tolerances (`AbsTol` and `RelTol`) that are applied to the integration-in-time algorithm performed by `ode15s`. For further details, please refer to the [Matlab documentation](https://mathworks.com/help/matlab/ref/odeset.html).
Note that, by default, the code follows a routine in which the solution procedure is attempted up to three times. On the second and third attempts (if needed), the code uses values for the two error tolerances that are 10 and 100 times smaller (or "stricter"), respectively.

If you encounter a different problem, please create an Issue on the GitHub website, add details of the problem (including the error message and MATLAB version number) and attach the parameters.m file in use when the problem occurred. For other enquiries, please contact N.E.Courtier(at)soton.ac.uk.


# How to Cite this Code

Please cite the release paper published in the [Journal of Computational Electronics](https://link.springer.com/article/10.1007/s10825-019-01396-2) by using the [citation.bib](citation.bib) file.
