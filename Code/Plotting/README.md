
This folder contains scripts for plotting the solution, after the solution structure `sol` is loaded from an output .mat file using a command such as:

`load('Data/simulation.mat');`

For example:
- plot_sections.m can be used to plot a J-V curve
- plot_recombination.m plots the rates of different types of bulk recombination
- plot_dstrbns.m generates plots of the dstrbns of P,n,p and phi as functions of space
- plot_IS.m can be used to plot data from any impedance simulation in Nyquist and Bode plots
- animate_sections.m renders animated plots of current, carrier distributions, electric potential and recombination rates during the simulation and saves the video as an mp4 file

- IonMonger_import.py is a Python file that imports a saved solution file into the Python environment and extracts most of the key variables for plotting or analysis