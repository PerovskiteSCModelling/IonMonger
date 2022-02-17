function params = parameters_template(iter)
% This function returns a parameter structure that contains all the inputs
% required for the simulation(s). In this function, the user should choose
% the values of the variables in the following sections:
% 1) Program settings
% 2) Parameter input
% 3) Simulation protocol
% The file should then be saved as parameters.m, before running master.m.
% To iterate over a number of different parameter sets, the user can make
% use of the optional integer input iter.

% Allow this function to be called without an argument
if nargin<1, iter=1; end

%% Program settings

% Options for the user
workfolder  = './Data/'; % the folder to which data will be saved,
% note that the string must end with a forward slash
OutputFcn   = 'PrintVolt'; % ode15s optional message function, choose
% either 'PrintVolt' or 'PrintTime' (which can be found in Code/Common/)
Verbose     = false; % set this to false to suppress message output,
% note that this option overwrites the previous two if Verbose=false
UseSplits   = true; % set this to false to make a single call to ode15s

% Resolution and error tolerances
N    = 400; % Number of subintervals, with N+1 being the number of grid points
rtol = 1e-6; % Relative temporal tolerance for ode15s solver
atol = 1e-8; % Absolute temporal tolerance for ode15s solver
phidisp = 100; % displacement factor for electric potential (VT) (optional)

%% Parameter input

% Physical constants
eps0  = 8.854187817e-12;  % permittivity of free space (Fm-1)
q     = 1.6021766209e-19; % charge on a proton (C)
Fph   = 1.4e21;           % incident photon flux (m-2s-1)
kB    = 8.61733035e-5;    % Boltzmann constant (eVK-1)

% Perovskite parameters
T     = 298;       % temperature (K)
b     = 400e-9;    % perovskite layer width (m) (normally between 150-600nm)
epsp  = 24.1*eps0; % permittivity of perovskite (Fm-1)
alpha = 1.3e7;     % perovskite absorption coefficient (m-1)
Ec    = -3.7;      % conduction band minimum (eV)
Ev    = -5.4;      % valence band maximum (eV)
Dn    = 1.7e-4;    % perovskite electron diffusion coefficient (m2s-1)
Dp    = 1.7e-4;    % perovskite hole diffusion coefficient (m2s-1)
gc    = 8.1e24;    % conduction band density of states (m-3)
gv    = 5.8e24;    % valence band density of states (m-3)

% Ion parameters
N0    = 1.6e25;        % typical density of ion vacancies (m-3)
Plim  = inf;           % limiting density of ion vacancies (m-3) (can choose inf)
D     = @(Dinf, EA) Dinf*exp(-EA/(kB*T)); % diffusivity relation
DIinf = 6.5e-8;        % high-temp. vacancy diffusion coefficient (m2s-1)
EAI   = 0.58;          % iodide vacancy activation energy (eV)
DI    = D(DIinf, EAI); % diffusion coefficient for iodide ions (m2s-1)

% Direction of light
inverted = false; % choose false for a standard architecture cell (light
% entering through the ETL), true for an inverted architecture cell
% (light entering through the HTL)

% ETL parameters
dE    = 1e24;    % effective doping density of ETL (m-3)
% EfE   = -4.1;    % doped electron quasi-Fermi level in ETL (eV) (optional, overrides dE)
gcE   = 5e25;    % effective conduction band DoS in ETL (m-3)
EcE   = -4.0;    % conduction band reference energy in ETL (eV)
bE    = 100e-9;  % width of ETL (m)
epsE  = 10*eps0; % permittivity of ETL (Fm-1)
% DE    = 1e-5;    % electron diffusion coefficient in ETL (m2s-1)
muE   = 3.8e-4;  % electron mobility in ETL (m2V-1s-1) (optional, overrides DE)
stats.ETL = struct('model','FermiDirac',... 
                    'Boltzmann',false); % ETL statistical model (choose model
% from 'FermiDirac', 'GaussFermi', or 'Blakemore')

% HTL parameters
% dH    = 1e24;    % effective doping density of HTL (m-3) 
EfH   = -4.852;    % doped hole quasi-Fermi in HTL (eV) (optional, overrides dH)
gvH   = iter;    % effective valence band DoS in HTL (m-3)
EvH   = -5.09;   % valence band reference energy in HTL (eV)
bH    = 200e-9;  % width of HTL (m)
epsH  = 3*eps0;  % permittivity of HTL (Fm-1)
% DH    = 1e-6;    % hole diffusion coefficient in HTL (m2s-1)
muH   = 3.8e-5;  % electron mobility in HTL (m2V-1s-1) (optional, overrides DH)
stats.HTL = struct('model','GaussFermi',... 
                    'Boltzmann',false,...
                    's',3.73); % HTL statistical model (choose model
% from 'FermiDirac', 'GaussFermi', or 'Blakemore')

% Metal contact parameters (optional)
% Ect   = -4.1;    % cathode workfunction (eV)
% Ean   = -5.0;    % anode workfunction (eV)

% Bulk recombination
tn    = 3e-9;    % electron pseudo-lifetime for SRH (s)
tp    = 3e-7;    % hole pseudo-lifetime for SRH (s)
beta  = 0e-10;       % bimolecular recombination rate (m3s-1)
Augn  = 0e-29;       % electron-dominated Auger recombination rate (m6s-1)
Augp  = 0e-29;       % hole-dominated Auger recombination rate (m6s-1)

% Interface recombination (max. velocity ~ 1e5)
betaE = 0e-22;       % ETL/perovskite bimolecular recombination rate (m3s-1)
betaH = 0e-22;       % perovskite/HTL bimolecular recombination rate (m3s-1)
vnE   = 1e5;     % electron recombination velocity for SRH (ms-1)
vpE   = 1e1;      % hole recombination velocity for SRH (ms-1)
vnH   = 1e-1;     % electron recombination velocity for SRH (ms-1)
vpH   = 1e5;     % hole recombination velocity for SRH (ms-1)

% Parasitic resistances (optional)
% Rs    = 0;       % external series resistance (Ohm)
% Rp    = Inf;     % parallel or shunt resistance (Ohm) (can choose Inf)
% Acell = 1;       % cell area (cm2) (only used to scale the series resistance)

%% Option to set initial distributions from a saved solution

% Name of a saved solution structure. Initial distributions will be taken
% from the final distributions of the saved solution. In this case,
% `applied_voltage` does not need an initial voltage.

% input_filename = 'SaveFiles/example1.mat';

%% Non-dimensionalise model parameters and save all inputs

% Compile all parameters into a convenient structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

% Non-dimensionalise the user-defined input parameters
params = nondimensionalise(params);

% Unpack variables and functions needed in the rest of this function
[tstar2t, psi2Vap, Upsilon, Vbi] = struct2array(params, ...
    {'tstar2t','psi2Vap','Upsilon','Vbi'});


%% Simulation protocol
% In order to make use of construct_protocol.m, instructions must be given
% in a specific order and format. Please see the GUIDE.md. Otherwise one
% can specify their own dimensionless functions of time (light and psi),
% dimensionless vectors (time and splits) and option whether to findVoc.

% Light protocol (either {a single value} or a protocol including an
% initial value, set to 1 for measurements in the light, 0 in the dark)
light_intensity = ...
    {1};

% Voltage protocol (either {'open-circuit'}, {a single value} or a protocol
% beginning with either 'open-circuit' or an initial value, in Volts). For
% impedance spectroscopy protocols, see GUIDE.md.
applied_voltage = ...
    {Vbi, ... % steady-state initial value
    'tanh', 5, 1.2, ...
    'linear', 1.2/1e-4, 0, ... % reverse scan
    'linear', 1.2/1e-4, 1.2 ... % reverse scan
    };

% Impedance protocol template:
% applied_voltage = {'impedance', ...
%     1e-4, ...   % minimum impedance frequency (Hz)
%     1e7, ...    % maximum impedance frequency (Hz)
%     0.9, ...  	% DC voltage (V)
%     20e-3, ...  % AC voltage amplitude (V)
%     200, ...      % number of frequencies to sample
%     5, ...      % number of sine waves
%     10};        % time to hold at DC voltage (s)

reduced_output = true; % set to true to reduce the amount of data retained ...
% in impedance simulations

% Choose whether the time points are spaced linearly or logarithmically
time_spacing = 'lin'; % set equal to either 'lin' (default) or 'log'

%% Create the simulation protocol and plot (if Verbose)

% Create the protocol and time points automatically, psi=(Vbi-Vap)/(2*kB*T)
[light, psi, time, splits, findVoc] = ...
    construct_protocol(params,light_intensity,applied_voltage,time_spacing);

% *** If defining one's own simulation protocol, define it here! ***

% Apply the options defined above
if inverted, inv = -1; else, inv = 1; end
if ~UseSplits, splits = time([1,end]); end

% Define the charge carrier generation function G(x,t)
G = @(x,t) light(t).*Upsilon./(1-exp(-Upsilon)).*exp(-Upsilon*(inv*x+(1-inv)/2));

% Plot the light regime
if Verbose
    if ishandle(98), clf(98); end; figure(98);
    plot(tstar2t(time),light(time),'Color',[0.93 0.69, 0.13]);
    xlabel('Time (s)'); ylabel('Light intensity (Sun equiv.)');
    title('light(t)');
    drawnow;
end

% Plot the voltage regime
if Verbose
    if ishandle(99), clf(99); end; figure(99);
    if isnan(psi(time(end)))
        title('Open-circuit');
    else
        plot(tstar2t(time),psi2Vap(psi(time)));
        xlabel('Time (s)'); ylabel('Applied Voltage (V)');
        title('V(t)');
        if findVoc
            hold on; plot(0,Vbi,'o','MarkerSize',8);
            title('V(t) except the voltage starts from Voc, not Vbi as shown here');
        end
        if strcmp(applied_voltage{1},'impedance')
            title({'example impedance protocol','at a frequency of 1Hz'})
            hold on
            plot(tstar2t(time(end-200:end)),psi2Vap(psi(time(end-200:end))),'r')
            legend({'','phase analysis region'},'Location','best')
            hold off
        end
        ylim([min([0,psi2Vap(psi(time))]), ceil(2*max(psi2Vap(psi(time))))/2])
    end
    drawnow;
end


%% Compile more parameters into the params structure
vars = setdiff(setdiff(who,fieldnames(params)),{'params','vars','i'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

% Make the folder in which to save the output data (specified above)
if exist(workfolder,'dir')~=7, mkdir(workfolder); end


end
