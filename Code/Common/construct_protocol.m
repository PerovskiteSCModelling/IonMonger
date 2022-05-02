function [light, psi, time, splits, findVoc] = ...
    construct_protocol(params,light_intensity,applied_voltage,time_spacing)
% This function can be used to construct the necessary voltage and light
% inputs from an applied_voltage protocol made of either a preconditioning
% voltage or the string 'open-circuit', optionally followed by a specific
% protocol described by a set of components consisting of:
%    {'type of slope', duration in seconds, voltage to end at}
%    [ Note that you can change scan rate to duration via:
%      Duration (s) = Voltage Change (V)/Scan Rate (V/s)  ]
% and a light_intensity protocol made of a preconditioning light intensity,
% optionally followed by an equivalent set of components:
%    {'type of slope', duration in seconds, intensity to end at}
% The choice of time_spacing can be either 'lin' for linearly-spaced time
% points or 'log' for logarithmically-spaced time points. This function is
% called by parameters.m. If both a voltage and an illumination protocol
% are specified, they must describe the same length of time. The outputs of
% the function are two functions of time (light and psi), two vectors of
% time points (time and splits) and the option whether to first findVoc.
% To create an impedance protocol, begin the applied_voltage cell with the
% string 'impedance', followed by the necessary protocol parameters (see
% user guide). In this case, the function will return a dummy protocol for
% a 1Hz impedance measurement for plotting purposes. 

% Parameter input
[Vbi, t2tstar, Vap2psi] = struct2array(params, {'Vbi','t2tstar','Vap2psi'});

if strcmp(applied_voltage{1},'impedance')
    % If the protocol is an impedance spectrum, create a dummy input for
    % plotting if Verbose
    
    V0 = applied_voltage{4}; % DC voltage
    Vp = applied_voltage{5}; % AC voltage amplitude
    t = 5;                  % time spent in steady state (s)
    n = applied_voltage{7};  % number of periods to simulate
    applied_voltage = {V0,'tanh',t,V0}; % protocol for steady state
    freq = 1; % example frequency (Hz)
    for i = 1:n
        % add n sin functions to the protocol
        applied_voltage{end+1} = 'sin';
        applied_voltage{end+1} = 1/freq;
        applied_voltage{end+1} = V0+Vp;
    end
end
    
% check for initial voltage
if length(applied_voltage)>1 & ~ischar(applied_voltage{2})
    % If second entry is not a character vector, the initial voltage has
    % been omitted
    if isfield(params,'input_filename')
        % check for specified initial conditions
        load(params.input_filename)
        
        % add the saved voltage to the beginning of the protocol
        applied_voltage = {sol.V(end), applied_voltage{:}}; 
    else
        error(['applied_voltage did not specify an initial voltage and '...
            'no specified initial distribution has been found.'])
    end
else
    if isfield(params,'input_filename')
        % initial voltage has been specified but an initial distribution
        % has also been specified
        warning(['Initial voltage was specified in applied_voltage but a ' ...
            'saved initial distribution has also been specified. This ' ...
            'value will override the initial voltage.'])
        load(params.input_filename)
        % replace initial voltage with saved voltage 
        applied_voltage{1} = sol.V(end); % replace initial voltage with saved voltage 
    end
end

if length(light_intensity)==1
    
    % If there is only a single input value for the light, set the light
    % intensity equal to this value for all time
    G_value = light_intensity{1};
    light = @(t) G_value*ones(size(t));
    G_splits = 0;
    
else
    
    % Find the strings in the cell array "light_intensity" to determine the
    % components of the protocol
    G_index = NaN(1,length(light_intensity));
    for i = 1:length(light_intensity)
        if ischar(light_intensity{i})
            G_index(i) = 1;
        end
    end
    G_index = find(G_index==1);
    
    % Define non-dimensional splits
    G_splits = t2tstar(cumsum([0, light_intensity{G_index+1}]));
    G_values = [light_intensity{[G_index(1:end)-1, end]}];
    G_compnts = {light_intensity{G_index}};

    % Define light as a function of time
    light = @(t) put_together(t, G_splits, G_values, G_compnts);
    
end

if length(applied_voltage)==1
    
    % If there is only a single input value for the voltage, it either
    % means simulating at open-circuit or a fixed voltage for all time
    if ischar(applied_voltage{1})
        psi = @(t) NaN; % value not used, just needs to be a function
        findVoc = true;
    else
        V_value = applied_voltage{1};
        psi = @(t) V_value*ones(size(t));
        findVoc = false;
    end
    V_splits = 0;
    
else
    
    % If the protocol starts from open-circuit, temporarily define the
    % first section of voltage protocol from the cell's built-in voltage,
    % rather than the as-yet-unknown Voc, and turn the findVoc option on
    if ischar(applied_voltage{1}) && strcmp(applied_voltage{1},'open-circuit')
        applied_voltage{1} = Vbi;
        findVoc = true;
    else
        findVoc = false;
    end
    
    % Find strings in the cell array "applied_voltage" to determine the components
    V_index = NaN(1,length(applied_voltage));
    for i = 1:length(applied_voltage)
        if ischar(applied_voltage{i})
            V_index(i) = 1;
        end
    end
    V_index = find(V_index==1);

    % Define non-dimensional splits
    V_splits = t2tstar(cumsum([0, applied_voltage{V_index+1}]));
    V_values = Vap2psi([applied_voltage{[V_index(1:end)-1, end]}]);
    V_compnts = {applied_voltage{V_index}};

    % Define psi as a function of time
    psi = @(t) put_together(t, V_splits, V_values, V_compnts);
    
end

% Merge the split times into one vector (which includes both the first
% and last time points)
if length(G_splits)==1 && length(V_splits)==1
    error('One of the protocols must include a length of time.');
elseif length(G_splits)>1 && length(V_splits)>1 && G_splits(end)~=V_splits(end)
    error('The protocols must describe the same length of time.');
else
    splits = unique([G_splits, V_splits]);
end

% Define time as a vector composed of the split times and 99 points evenly
% distributed between one split time and the next
time = NaN(1,(length(splits)-1)*100+1);
for i = 1:length(splits)-1
    if strcmp(time_spacing,'lin') % For linearly spaced points in time
        time((i-1)*100+1:i*100+1) = linspace(splits(i),splits(i+1),101);
    elseif strcmp(time_spacing,'log') % For logarithmic spacing
        if i==1
            time(1) = 0; % must treat t=0 separately
            time(2:101) = logspace(log10(splits(2)/100),log10(splits(2)),100);
        else
            time((i-1)*100+1:i*100+1) = logspace(log10(splits(i)),log10(splits(i+1)),101);
        end
    else
        error('TimeSpacing must equal one of the two options ''lin'' or ''log''.');
    end
end

end

function out = put_together(t, splits, values, compnts)
% Constructs a function of time from the set of instructions provided
% Preallocate
out = NaN(size(t));
% Define the value of the function by each component
for i = 1:length(splits)-1
    eval([ ...
        'idxs = splits(i)<=t & t<splits(i+1);', ...
        sprintf('out(idxs) = part_%s',compnts{i}),...
        '(t(idxs), splits(i), splits(i+1), values(i), values(i+1));' ...
        ]);
    if strcmp(compnts{i},'sin')
        values(i+1) = values(i);
    end
end
% Afterwards, remain at the final value
idxs = splits(end)<=t;
out(idxs) = values(end);
end


%% Functions used by the above code:

function out = part_linear(t, t_start, t_end, V_start, V_end)
% Linear transition from [t_start,V_start] to [t_end,V_end]
out = V_start + (V_end-V_start) .* (t-t_start)./(t_end-t_start);
end

function out = part_cosine(t, t_start, t_end, V_start, V_end)
% Uses the cosine function to transition from V_start to V_end in a
% sigmoidal fashion between t_start and t_end
out = V_start+(V_end-V_start).* ...
    (0.5-0.5.*cos(pi*(t-t_start)./(t_end-t_start)));
end

function out = part_sin(t, t_start, t_end, V_start, V_max)
% Uses the sin function between V_max and V_start-(V_max-V_start) in a
% single sinusoid between t_start and t_end.
% Note that V_end has been replaced by V_max.
out = V_start+(V_max-V_start).* ...
    sin(2*pi*(t-t_start)/(t_end-t_start));
end

function out = part_tanh(t, t_start, t_end, V_start, V_end)
% Uses the tanh function to transition from V_start to V_end via a sharp
% rise which then flattens off between t_start and t_end, with steepness k
k = 3;
out = V_start+(V_end-V_start).* ...
    tanh(k*(t-t_start)./(t_end-t_start))/tanh(k);
end
