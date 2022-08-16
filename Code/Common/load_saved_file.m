function sol_start = load_saved_file(params,vectors)
% This function creates a vector containing a consistent set of initial
% conditions by interpolating an existing solution onto the current
% spatial grid. The inputs are structures containing the parameters and
% vectors needed to nondimensionalise and interpolate the solution.

% Unpack parameters and vectors
[b, N0, VT, dE, dH, n0, p0, Verbose] = ...
    struct2array(params,{'b','N0','VT','dE','dH','n0','p0','Verbose'});
[x, xE, xH] = struct2array(vectors,{'x','xE','xH'});
b = b*1e9; % Change the width b in metres to nanometres

if Verbose
    disp(['Using initial distributions from the saved file ' ... 
        params.input_filename])
end

% Load file
load(params.input_filename);

% Unpack dimensional spatial vectors from input file
[old_x, old_xE, old_xH] = struct2array(sol.vectors,{'x','xE','xH'});

% Set the chosen time point as the final time by default
tt = size(sol.dstrbns.phi,1);

% Unpack the solution at the chosen time point
P    = sol.dstrbns.P(tt,:);
phi  = sol.dstrbns.phi(tt,:);
n    = sol.dstrbns.n(tt,:);
p    = sol.dstrbns.p(tt,:);
phiE = sol.dstrbns.phiE(tt,:);
nE   = sol.dstrbns.nE(tt,:);
phiH = sol.dstrbns.phiH(tt,:);
pH   = sol.dstrbns.pH(tt,:);

% Interpolate onto new spatial grid and nondimensionalise
P    = interp1(old_x, P,   x*b )/N0;
phi  = interp1(old_x, phi, x*b )/VT;
n    = interp1(old_x, n,   x*b )/n0;
p    = interp1(old_x, p,   x*b )/p0;
phiE = interp1(old_xE,phiE,xE*b)/VT;
nE   = interp1(old_xE,nE,  xE*b)/dE;
phiH = interp1(old_xH,phiH,xH*b)/VT;
pH   = interp1(old_xH,pH,  xH*b)/dH;

% Eliminate superfluous points for the electric potential at the interfaces
phiE = phiE(1:end-1);
phiH = phiH(2:end);

% Combine the distributions into one vector
sol_start = [P; phi; n; p; phiE; nE; phiH; pH];

if any(isnan(sol_start))
    error(['There was an error in interpolating the saved ' ... 
        'distributions onto the new spatial grid.']);
end
    
end
