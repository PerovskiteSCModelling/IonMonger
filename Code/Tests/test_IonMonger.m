function tests = test_IonMonger

    % Add IonMonger directories to path
    cd('../../');
    addpath(cd);
    reset_path();
    cd('./Code/Tests/');

    % Find and run tests
    tests = functiontests(localfunctions);
    
end

function regression_test(testCase)
% Test to check whether the functionality of IonMonger v1.0 is retained in
% the current version, by comparing the output for the example simulation
% presented in Figure 1(a) of the release paper, which can be found at:
% https://doi.org/10.1007/s10825-019-01396-2.

% Load the original data and parameters structure
load('./Data/simulation_1.0.mat','sol');
params = sol.params;

% Update and import parameters
params.Verbose = false;
VT = params.VT;
jay = params.jay;

% Update parameters structure with new parameters
params = nondimensionalise(params);

% Solve the equations
newsol = numericalsolver(params);

% Define absolute and relative difference functions
abs_diff = @(new,old) max(abs(new-old));
rel_diff = @(new,old) max(max(abs(new-old)./max(abs(old),abs(mean(old,2)))));

% Compare the results
P_diff    = rel_diff(newsol.dstrbns.P   ,sol.dstrbns.P   );
phi_diff  = abs_diff(newsol.dstrbns.phi ,sol.dstrbns.phi )/VT;
n_diff    = rel_diff(newsol.dstrbns.n   ,sol.dstrbns.n   );
p_diff    = rel_diff(newsol.dstrbns.p   ,sol.dstrbns.p   );
phiE_diff = abs_diff(newsol.dstrbns.phiE,sol.dstrbns.phiE)/VT;
nE_diff   = rel_diff(newsol.dstrbns.nE  ,sol.dstrbns.nE  );
phiH_diff = abs_diff(newsol.dstrbns.phiH,sol.dstrbns.phiH)/VT;
pH_diff   = rel_diff(newsol.dstrbns.pH  ,sol.dstrbns.pH  );
V_diff    = abs_diff(newsol.V ,sol.V )/max(VT ,abs(sol.V ));
J_diff    = abs_diff(newsol.J ,sol.J )/max(jay,abs(sol.J ));
Jl_diff   = abs_diff(newsol.Jl,sol.Jl)/max(jay,abs(sol.Jl));
Jr_diff   = abs_diff(newsol.Jr,sol.Jr)/max(jay,abs(sol.Jr));
max_diff  = max([P_diff, phi_diff, n_diff, p_diff, ...
                 phiE_diff, nE_diff, phiH_diff, pH_diff, ...
                 V_diff, J_diff, Jl_diff, Jr_diff]);

% Check the difference is less than the tolerance
verifyLessThan(testCase, max_diff, 5e-5);

end

function integration_test(testCase)
% Test to check whether any routine errors are returned when simulating
% different combinations of light and voltage protocols for the default
% parameter set. This test does not check the output.

% Create a structure filled with the default parameters
params = parameters_template();
params.Verbose = false;
    
% Set standard resolution and error tolerances
params.N = 400;
params.rtol = 1e-6;
params.atol = 1e-10;

% Define the experimental protocols
light_protocol = {};
voltage_protocol = {};
% 1. Example J-V scan from Figure 1(a) of the release paper in the light
light_protocol{end+1} = {1};
voltage_protocol{end+1} = {1.2, 'linear', 1.2, 0, 'linear', 1.2, 1.2};
% 2. Switching the light off and on while at open-circuit
light_protocol{end+1} = {1, 'linear', 1, 0, 'linear', 1, 1};
voltage_protocol{end+1} = {'open-circuit'};

% Test each protocol in turn
error_count = 0;
for i = 1:length(light_protocol)
    try
        
        % Select the experimental protocol
        light_intensity = light_protocol{i};
        applied_voltage = voltage_protocol{i};
        time_spacing = 'lin'; % set equal to either 'lin' (default) or 'log'
        
        % Create the protocol, time points and generation function
        [params.light, params.psi, params.time, params.splits, params.findVoc] = ...
            construct_protocol(params,light_intensity,applied_voltage,time_spacing);
        [light, inv, Upsilon]= struct2array(params,{'light','inv','Upsilon'});
        params.G = @(x,t) light(t).*Upsilon./(1-exp(-Upsilon)).*exp(-Upsilon*(inv*x+(1-inv)/2));
        
        % Solve the equations
        numericalsolver(params);
        
    catch
        error_count = error_count + 1;
    end
    
end

% Check that there are no errors
verifyEqual(testCase, error_count, 0);

end
