function [J, Jl, Jr] = calculate_currents(params,vectors,dstrbns)
% This function computes the values of the total current density J, the
% current density lost to recombination across the ETL/perovskite interface
% Jl, and the current density lost to recombination across the
% perovskite/HTL interface Jr. The inputs are structures containing the
% input paramters, spatial vectors and solution variables. The outputs are
% three column vectors containing the dimensionless values of J, Jl and Jr
% at each point in time.

% Unpack relevant parameters, vectors and solution variables
[N, time, Kn, Kp, dpt, dpf, Rl, Rr, ARp] ...
    = struct2array(params,{'N','time','Kn','Kp','dpt','dpf','Rl','Rr','ARp'});
dx  = vectors.dx;
[P, phiE, phi, phiH, n, p] ...
    = struct2array(dstrbns, {'P','phiE','phi','phiH','n','p'});

% Define necessary vectors and the indice closest to the mid-point of the
% perovskite layer
dt  = diff(time)';
TT  = (1:length(time))';
mid = ceil((N+1)/2);

% Define the dimensionless charge carrier and displacement current densities
% at t=T and the half-point x(k)+dx(k)/2 (i.e. at the half-point denoted by j=k)
jn = @(T,k)  Kn./dx(k).*(n(T,k+1)-n(T,k)-(n(T,k+1)+n(T,k)).*(phi(T,k+1)-phi(T,k))./2);
jp = @(T,k) -Kp./dx(k).*(p(T,k+1)-p(T,k)+(p(T,k+1)+p(T,k)).*(phi(T,k+1)-phi(T,k))./2);
jd = @(T,k)  dpt./dx(k).*(phi(T,k+1)-phi(T,k)-phi(T-1,k+1)+phi(T-1,k))./dt(T-1);
jf = @(T,k) -dpf./dx(k).*(P(T,k+1)-P(T,k)+(P(T,k+1)+P(T,k)).*(phi(T,k+1)-phi(T,k))./2);
    
% Evaluate the dimensionless current densities at the mid-point of the perovskite
Jn = jn(TT,mid); % dimensionless electron current density
Jp = jp(TT,mid); % dimensionless hole current density
Jd = [NaN; jd(TT(2:end),mid)]; % dimensionless displacement current density
Jf = jf(TT,mid); % dimensionless ionic flux displacement current density
Js = (params.pbi-(phiE(TT,1)-phiH(TT,end)))/ARp; % loss due to shunt resistance

% Calculate the total dimensionless photocurrent density
J = Jn+Jp-Jd-Jf-Js;

% Calculate dimensionless losses from interfacial recombination
Jl = -Rl(n(TT,1),p(TT,1)); % (negative) current density loss at ETL interface
Jr = -Rr(n(TT,N+1),p(TT,N+1)); % (negative) current density loss at HTL interface

end
