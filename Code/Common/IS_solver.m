function sol = IS_solver(base_params)

nf = base_params.applied_voltage{7}; % number of frequencies to be sampled
min_f = base_params.applied_voltage{2};
max_f = base_params.applied_voltage{3};
n_wave = base_params.applied_voltage{8};

freqs = logspace(log10(min_f),log10(max_f),nf);

% find steady state at the DC voltage
fprintf('solving for steady state conditions at DC voltage \n')
params = base_params;
V0 = base_params.applied_voltage{4};
t = base_params.applied_voltage{6};
params.applied_voltage = {params.Vbi,'tanh',t,V0};

[params.light, params.psi, params.time, params.splits, params.findVoc] = ...
    construct_protocol(params,params.light_intensity,params.applied_voltage,params.time_spacing);

sol_init = numericalsolver(params);
savestr = 'Data/DC_sol';
save_end_state(sol_init,savestr)

parfor j = 1:nf
    disp(['solving frequency ' num2str(j) '/' num2str(nf)])
    params = base_params;
    params.Verbose = false;
    params.UseSplits = false;
    params.input_filename = savestr;
    
    V0 = base_params.applied_voltage{4};
    Vp = base_params.applied_voltage{5};
    t = base_params.applied_voltage{6};
    
    % construct initial stabilisation
    params.applied_voltage = {'linear',t,V0};
    for k = 1:n_wave % add multiple waves
        params.applied_voltage{end+1} = 'sin';
        params.applied_voltage{end+1} = 1/freqs(j);
        params.applied_voltage{end+1} = V0+Vp;
    end
    
    [params.light, params.psi, params.time, params.splits, params.findVoc] = ...
    construct_protocol(params,params.light_intensity,params.applied_voltage,params.time_spacing);

    try
        sol = numericalsolver(params);
        sols(j) = sol;
    catch me
        warning(['frequency ' num2str(j) ' encountered an error'])
        disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end
end

for j = 1:length(sols)
    sols(j).impedance_protocol = base_params.applied_voltage; % retain overall protocol
end

% analyse output
[X,R] = impedance_analysis(sols);

%% decide what information to retain in the impedance sol structure

if base_params.reduced_output
    sol.params = base_params;
    sol.X = X;
    sol.R = R;
    sol.freqs = freqs;
else
    sol = sols;
end

end
