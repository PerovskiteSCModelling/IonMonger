function sol_init = initial_conditions(psi0,params,vectors,matrices)
% This function creates a vector containing a steady-state solution to the 
% initial problem. The inputs are structures containing the parameters,
% vectors and matrices needed by the solver.

% Parameter input
[chi, nc, pc, Verbose, N, NE, NH, phidisp] = struct2array(params, {'chi','nc', ...
    'pc','Verbose', 'N', 'NE', 'NH', 'phidisp'});
[x, xE, xH] = struct2array(vectors, {'x','xE','xH'});

% Define uniform profiles for the ion vacancy density and electric potential
P_init    = ones(size(x));
phi_init  = zeros(size(x))+phidisp;
phiE_init = zeros(size(xE(1:NE)))+phidisp;
phiH_init = zeros(size(xH(1:NH)))+phidisp;

% Compute profiles for the carrier concentrations from a quasi-steady BVP
y_guess = bvpinit(x',@(x) [x+(1-x)/chi; 0*x; chi*x+(1-x); 0*x]);
sol = bvp4c(@(x,y) yode(x,y,params),@(ya,yb) ybcs(ya,yb,params),y_guess);
solx = deval(sol,x'); p_init = solx(1,:)'; n_init = solx(3,:)';

% Define tanh profiles for the carrier concentrations across the TLs
stE = 3/xE(1);
nE_init = nc+(1-nc)*tanh(stE*(xE(1)-xE))/tanh(stE*(xE(1)-xE(end)));
stH = 3/(xH(end)-1);
pH_init = pc+(1-pc)*tanh(stH*(xH(end)-xH))/tanh(stH*(xH(end)-xH(1)));

% Combine the initial conditions into one vector to pass to fsolve
sol_init  = [P_init; phi_init; n_init; p_init; ... % perovskite
             phiE_init; nE_init; ... % electron transport layer
             phiH_init; pH_init]; % hole transport layer
          
% Define the settings for the call to fsolve
fsoptions = optimoptions('fsolve','MaxIterations',20);
if Verbose, fsoptions.Display = 'iter'; else, fsoptions.Display = 'off'; end

% Use the initial guess to obtain an approximate steady-state solution
if exist('AnJac.m','file')
    fsoptions.SpecifyObjectiveGradient = true;
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS_AnJac(u,psi0, ...
        params,vectors,matrices,'init'),sol_init,fsoptions);
else
    fsoptions.JacobPattern = Jac(params,'init');
    [sol_init,~,exitflag,~] = fsolve(@(u) RHS(0,u,psi0, ...
        params,vectors,matrices,'init'),sol_init,fsoptions);
end
if exitflag<1
    warning(['Steady-state initial conditions could not be found to ' ...
        'a high degree of accuracy and may be unphysical.']);
end

% Ensure all the algebraic equations are satisfied as exactly as possible
sol_init = apply_Poisson(sol_init,params,vectors,matrices);

end


% Quasi-steady BVP for the carrier concentrations
% y(1) = p(x), y(2) = jp(x) = -Kp*dp/dx, y(3) = n(x), y(4) = jn(x) = Kn*dn/dx
function dpdx = yode(x,y,params)
[G, R, Kp, Kn] = struct2array(params,{'G','R','Kp','Kn'});
dpdx = [-y(2)/Kp; ...
        G(x,0)-R(y(3),y(1),1); ...
        y(4)/Kn; ...
        -(G(x,0)-R(y(3),y(1),1))];
end
function res = ybcs(ya,yb,params)
[Rl, Rr, AH, AE, omegaH, omegaE, SHinv, SEinv] = struct2array(params,{...
    'Rl','Rr', 'AH', 'AE', 'omegaH', 'omegaE', 'SHinv', 'SEinv'});
res = [yb(1)-AH/omegaH*exp(SHinv(omegaH)); ...
       ya(2)+Rl(1,ya(1)); ...
       ya(3)-AE/omegaE*exp(SEinv(omegaE)); ...
       yb(4)+Rr(yb(3),1)];
end
