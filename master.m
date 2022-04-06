% This is the master script for running a single simulation.

% Begin
clear;
tic;
fprintf('Computation started at %s\n', datestr(now));
reset_path();

% Create a structure filled with the user-defined parameters

params = parameters_template();

% Solve the equations
if strcmp(params.applied_voltage{1},'impedance') % check if voltage protocol ...
% specifies an impedance protocol
    sol = IS_solver(params);
else
    sol = numericalsolver(params);
end

% Save the data
save([params.workfolder,'simulation.mat'],'sol');

% Stop stopwatch and output nice message
fprintf('Completed simulation at %s, taking %s\n', ...
    datestr(now), secs2hms(toc) );

