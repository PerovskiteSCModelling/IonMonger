function params = nondimensionalise(params)
% This function automatically non-dimensionalises all the model parameters
% as required. The input and output both take the form of a structure.

% Parameter input
[N, q, Fph, kB, T, b, epsp, alpha, Ec, Ev, Dn, Dp, gc, gv, N0, DI, EcE, dE, ...
    gcE, bE, epsE, DE, EvH, dH, gvH, bH, epsH, DH, tn, tp, beta, Augn, Augp, ...
    betaE, betaH, vnE, vpE, vnH, vpH, Ect, Ean, Rs, Rp, Acell] ...
    = struct2array(params, ...
    {'N','q','Fph','kB','T','b','epsp','alpha','Ec','Ev','Dn','Dp','gc', ...
    'gv','N0','DI','EcE','dE','gcE','bE','epsE','DE','EvH','dH','gvH','bH', ...
    'epsH','DH','tn','tp','beta','Augn','Augp','betaE','betaH','vnE','vpE', ...
    'vnH','vpH','Ect','Ean','Rs','Rp','Acell'});

% Energy level parameters
VT = kB*T; % thermal voltage (V)
EfE = EcE-VT*log(gcE/dE); % workfunction of ETL (eV)
EfH = EvH+VT*log(gvH/dH); % workfunction of HTL (eV)
if ~any(Ect), Ect = EfE; end % cathode workfunction (eV)
if ~any(Ean), Ean = EfH; end % anode workfunction (eV)
Vbi = Ect-Ean; % built-in voltage (V)
pbi = Vbi/VT;  % non-dim. built-in voltage

% Perovskite parameters
Eg      = Ec-Ev;                % bandgap (eV)
LD      = sqrt(VT*epsp/(q*N0)); % Debye length (m)
lambda  = LD/b;                 % Debye length parameter
lam2    = lambda^2;             % Debye length parameter squared
ni      = sqrt(gc*gv)*exp(-Eg/(2*VT)); % intrinsic carrier density (m-3)
n0      = gc*exp((EfE-Ec)/VT);  % typical electron density in perovskite (m-3)
p0      = gv*exp((Ev-EfH)/VT);  % typical hole density in perovskite (m-3)
delta   = n0/N0;          % ratio of typical electron and ion densities
chi     = p0/n0;          % ratio of typical hole and electron densities
G0      = (Fph./b).*(1-exp(-alpha*b)); % typical rate of photogeneration (m-3s-1)
if nnz(DI) % the ion diffusion coefficient is non-zero
    Tion = b/DI*sqrt(VT*epsp/(q*N0));  % characteristic ionic timescale (s)
else
    Tion = n0/G0;                      % characteristic electronic timescale (s)
end
sigma   = n0/(G0*Tion);   % ratio of carrier and ionic timescales
Kp      = Dp*p0/(G0*b^2); % hole current parameter
Kn      = Dn*n0/(G0*b^2); % electron current parameter
Upsilon = alpha*b;        % parameter for Beer-Lambert Law photo-generation
jay     = q*G0*b/10; % rescaling factor for current density to be in mAcm-2
dpt     = epsp*VT/(q*G0*b^2*Tion); % displacement current density factor
dpf     = DI*N0/(G0*b^2); % ionic flux displacement current density factor

% Transport layer parameters
wE = bE/b;             % relative width of ETL
wH = bH/b;             % relative width of HTL
KE = DE*dE/(G0*b^2);   % ETL electron current parameter
KH = DH*dH/(G0*b^2);   % HTL hole current parameter
rE = epsE/epsp;        % relative ETL permittivity
rH = epsH/epsp;        % relative HTL permittivity
lamE2 = rE*N0/dE*lam2; % relative ETL Debye length parameter squared
lamE  = sqrt(lamE2);   % relative ETL Debye length parameter
lamH2 = rH*N0/dH*lam2; % relative HTL Debye length parameter squared
lamH  = sqrt(lamH2);   % relative HTL Debye length parameter
OmegaE = sqrt(N0/(rE*dE)); % ETL charge density parameter
OmegaH = sqrt(N0/(rH*dH)); % HTL charge density parameter

% Interface parameters
kE = n0/dE; % ratio between electron densities across ETL/perovskite interface
kH = p0/dH; % ratio between hole densities across perovskite/HTL interface

% Contact parameters
nc = gcE*exp((Ect-EcE)/VT)/dE; % non-dim. electron density at cathode interface
pc = gvH*exp((EvH-Ean)/VT)/dH; % non-dim. hole density at anode interface

% Check for Auger recombination parameters
if isempty(Augn) || isempty(Augp)
    [Augn, Augp] = deal(0); % no Auger recombination
    warning(['The parameters for Auger recombination (Augn and Augp) have not' ...
        ' been specified in the parameters file, so they have been set to zero.']);
end

% Bulk recombination parameters
ni2   = ni^2/(n0*p0);  % non-dim. n_i^2
brate = beta*n0*p0/G0; % rate constant for bimolecular recombination
Cn = Augn*n0^2*p0/G0; % Auger recombination coefficient
Cp = Augp*n0*p0^2/G0; % Auger recombination coefficient
if tp>0 && tn>0
    gamma   = p0/(tp*G0); % rate constant for SRH recombination
    tor     = tn*p0/(tp*n0); % ratio of SRH carrier lifetimes
    tor3    = (tn+tp)*ni/(tp*n0); % constant from deep trap approximation
else
    [gamma, tor, tor3] = deal(0); % no bulk SRH
end
% Auger recombination rate
Auger = @(n,p,Cn,Cp,ni2) (Cn*n+Cp*p).*(n.*p-ni2);
% SRH recombination rate (written in a way that should reduce numerical inaccuracy)
SRH = @(n,p,gamma,ni2,tor,tor3) ...
        gamma*(p-ni2./n)./(1+tor*p./n+tor3./n).*(n>=tor*p).*(n>=tor3) ...
      + gamma*(n-ni2./p)./(n./p+tor+tor3./p).*(tor*p>n).*(tor*p>tor3) ...
      + gamma*(p.*n-ni2)./(n+tor*p+tor3).*(tor3>n).*(tor3>=tor*p);
% Total bulk recombination rate
R = @(n,p,P) brate*(n.*p-ni2) ... % bimolecular
           + Auger(n,p,Cn,Cp,ni2) ... % Auger recombination
           + SRH(n,p,gamma,ni2,tor,tor3); % SRH recombination via trap states

% Interface recombination parameters
brateE = betaE*dE*p0/(b*G0); % rate constant for bimolecular recombination
brateH = betaH*n0*dH/(b*G0); % rate constant for bimolecular recombination
if vpE>0 && vnE>0
    gammaE = p0*vpE/(b*G0); % rate constant for SRH recombination
    torE   = p0*vpE/(dE*vnE); % ratio of SRH carrier lifetimes
    torE3  = (1/n0+vpE/(dE*vnE))*ni; % constant from deep trap state approximation
else
    [gammaE, torE, torE3] = deal(0); % no ETL/perovskite interface recombination
end
if vnH>0 && vpH>0
    gammaH = n0*vnH/(b*G0); % rate constant for SRH recombination
    torH   = n0*vnH/(dH*vpH); % ratio of SRH carrier lifetimes
    torH3  = (1/p0+vnH/(dH*vpH))*ni; % constant from deep trap state approximation
else
    [gammaH, torH, torH3] = deal(0); % no perovskite/HTL interface recombination
end
% Total ETL/perovskite interface recombination rate
Rl = @(n,p) brateE*(n.*p-ni2) ... % bimolecular
          + SRH(n,p,gammaE,ni2,torE,torE3); % SRH
% Total perovskite/HTL interface recombination rate
Rr = @(n,p) brateH*(n.*p-ni2) ... % bimolecular
          + SRH(p,n,gammaH,ni2,torH,torH3); % SRH

% Spatial grid parameters (consistent with the choice of N above)
X = 0.2; % percentage of grid points within an ionic Debye length of the interface
tanhfun = @(x,st) (tanh(st*(2*x-1))/tanh(st)+1)/2;
options = optimoptions(@fsolve,'Display','Off','OptimalityTolerance',1e-8);
st = fsolve(@(st) lambda-tanhfun(X,st),2,options);
A = @(b) (tanh(st*(1-2/N))-(1-double(b))*tanh(st))/double(b);
NE = round(2/(1-atanh(A(wE))/st)); % number of subintervals in the ETL
NH = round(2/(1-atanh(A(wH))/st)); % number of subintervals in the HTL

% Little functions to de/re-dimensionalise
TkT     = 2*kB*T; % the factor of 2 is because +-psi/2 at left/right
t2tstar = @(t) t./Tion;
tstar2t = @(tstar) tstar.*Tion;
Vap2psi = @(Vap) (Vbi-Vap)./TkT;
psi2Vap = @(psi) Vbi-psi.*TkT;

% External parameters
if isempty(Rs), Rs = 0; % default is zero series resistance
    disp('Assumming there is no series resistance.'); end
if ~any(Rp), Rp = Inf;  % default is infinite shunt resistance
    disp('Assumming there is infinite shunt resistance.'); end
if ~any(Acell), Acell = 1; end % default is cell area of 1 cm2
ARs = Rs*Acell/1e4*q*G0*b/VT; % non-dim. external series resistance x cell area
ARp = Rp*Acell/1e4*q*G0*b/VT; % non-dim. parallel/shunt resistance x cell area
Rsp = Rs/Rp; % ratio between series and parallel/shunt resistance

% Compile all parameters into the params structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

end
