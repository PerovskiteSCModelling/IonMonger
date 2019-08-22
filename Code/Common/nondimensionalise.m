function params = nondimensionalise(params)
% This function automatically non-dimensionalises all the model parameters
% as required. The input and output both take the form of a structure.

% Parameter input
[N, q, Fph, kB, T, b, epsp, alpha, Ec, Ev, Dn, Dp, gc, gv, N0, DI, EcE, dE, ...
    gcE, bE, epsE, DE, EvH, dH, gvH, bH, epsH, DH, tn, tp, beta, betaE, betaH, ...
    vnE, vpE, vnH, vpH] = struct2array(params, ...
    {'N','q','Fph','kB','T','b','epsp','alpha','Ec','Ev','Dn','Dp','gc', ...
    'gv','N0','DI','EcE','dE','gcE','bE','epsE','DE','EvH','dH','gvH','bH', ...
    'epsH','DH','tn','tp','beta','betaE','betaH','vnE','vpE','vnH','vpH'});

% Energy level parameters
VT = kB*T; % thermal voltage (V)
EfE = EcE-VT*log(gcE/dE); % workfunction of ETL (eV)
EfH = EvH+VT*log(gvH/dH); % workfunction of HTL (eV)
Vbi = EfE-EfH; % built-in voltage (V)
pbi = Vbi/VT;  % non-dim. built-in voltage

% Perovskite parameters
Eg      = Ec-Ev;                % bandgap (eV)
LD      = sqrt(VT*epsp/(q*N0)); % Debye length (m)
lambda  = LD/b;                 % Debye length parameter
lam2    = lambda^2;             % Debye length parameter squared
ni      = sqrt(gc*gv)*exp(-Eg/(2*VT)); % intrinsic carrier density (m-3)
delta   = dE/N0;          % ratio of typical electron and ion densities
chi     = dH/dE;          % ratio of typical hole and electron densities
Tion    = b/DI*sqrt(VT*epsp/(q*N0));   % characteristic ionic timescale (s)
G0      = (Fph./b).*(1-exp(-alpha*b)); % typical rate of photogeneration (m-3s-1)
sigma   = dE/(G0*Tion);   % ratio of carrier and ionic timescales
Kp      = Dp*dH/(G0*b^2); % hole current parameter
Kn      = Dn*dE/(G0*b^2); % electron current parameter
Upsilon = alpha*b;        % parameter for Beer-Lambert Law photo-generation
jay     = q*G0*b/10; % rescaling factor for current density to be in mAcm-2
dpt     = epsp*VT/(q*G0*b^2*Tion); % displacement current density factor
dpf     = DI*N0/(G0*b^2); % ionic flux displacement current denisty factor

% Transport layer parameters
wE = bE/b;      % relative width of ETL
wH = bH/b;      % relative width of HTL
KE = DE*Kn/Dn;  % ETL electron current parameter
KH = DH*Kp/Dp;  % HTL hole current parameter
rE = epsE/epsp; % relative ETL permittivity
rH = epsH/epsp; % relative HTL permittivity
lamE2 = rE*N0/dE*lam2; % relative ETL Debye length parameter squared
lamE  = sqrt(lamE2);   % relative ETL Debye length parameter
lamH2 = rH*N0/dH*lam2; % relative HTL Debye length parameter squared
lamH  = sqrt(lamH2);   % relative HTL Debye length parameter
OmegaE = sqrt(N0/(rE*dE)); % ETL charge density parameter
OmegaH = sqrt(N0/(rH*dH)); % HTL charge density parameter

% Interface parameters
kE = gc/gcE*exp((EcE-Ec)/VT); % ratio between electron densities across ETL/perovskite interface
kH = gv/gvH*exp((Ev-EvH)/VT); % ratio between hole densities across perovskite/HTL interface
n0 = kE*dE; % typical electron density in perovskite (m-3)
p0 = kH*dH; % typical hole density in perovskite (m-3)

% Bulk recombination parameters
ni2   = ni^2/(dE*dH);  % non-dim. n_i^2
brate = beta*dE*dH/G0; % rate constant for bimolecular recombination
if tp>0 && tn>0
    gamma   = dH/(tp*G0); % rate constant for hole-dominated recombination
    tor     = tn*dH/(tp*dE); % ratio of SRH carrier lifetimes
    tor3    = (tn+tp)*ni/(tp*dE); % constant from deep trap approximation
else
    [gamma, tor, tor3] = deal(0); % no bulk recombination
end
% Bulk recombination rate
R = @(n,p,AP) brate*(n.*p-ni2) ... % bimolecular
            + gamma*(p.*n-ni2)./(n+tor*p+tor3); % SRH

% Interface recombination parameters
brateE = betaE*dE*dH/(b*G0); % rate constant for bimolecular recombination
brateH = betaH*dE*dH/(b*G0); % rate constant for bimolecular recombination
if vpE>0 && vnE>0
    gammaE = dH*vpE/(b*G0); % rate constant for hole-dominated recombination
    torE   = dH*vpE/(dE*vnE); % ratio of SRH carrier lifetimes
    torE3  = (1/kE+vpE/vnE)*ni/dE; % constant from deep trap state approximation
else
    [gammaE, torE, torE3] = deal(0); % no ETL/perovskite interface recombination
end
if vnH>0 && vpH>0
    gammaH = dE*vnH/(b*G0); % rate constant for hole-dominated recombination
    torH   = dE*vnH/(dH*vpH); % ratio of SRH carrier lifetimes
    torH3  = (1/kH+vnH/vpH)*ni/dH; % constant from deep trap state approximation
else
    [gammaH, torH, torH3] = deal(0); % no perovskite/HTL interface recombination
end
% ETL/perovskite interface recombination rate
Rl = @(nE,p) brateE*(nE.*p-ni2/kE) ... % bimolecular
           + gammaE*(p.*nE-ni2/kE)./(nE+torE*p+torE3); % SRH
% perovskite/HTL interface recombination rate
Rr = @(n,pH) brateH*(n.*pH-ni2/kH) ... % bimolecular
           + gammaH*(pH.*n-ni2/kH)./(pH+torH*n+torH3); % SRH

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

% Compile all parameters into the params structure
vars = setdiff(who,{'params','vars'});
for i=1:length(vars), params.(vars{i}) = eval(vars{i}); end

end
