% This is the master script for running a single simulation.

% Begin
clear;
tic;
fprintf('Computation started at %s\n', datestr(now));
reset_path();

% Create a structure filled with the user-defined parameters
params = parameters();

% Solve the equations
sol = numericalsolver(params);

% Save the data
save([params.workfolder,'simulation.mat'],'sol');

% Stop stopwatch and output nice message
fprintf('Completed simulation at %s, taking %s\n', ...
    datestr(now), secs2hms(toc) );
