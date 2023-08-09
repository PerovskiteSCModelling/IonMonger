function [Q, Vs] = ionic_charge(params)
% A function to compute the Debye layer charge density over time.

% Unpack parameters
[psi, time, OmegaE, OmegaH, q, N0, epsp, VT, Vdc, Vap2psi] = ...
    struct2array(params,{'psi','time','OmegaE','OmegaH','q','N0', ...
                         'epsp','VT','Vdc','Vap2psi'});

% Solver settings
optimopts = optimoptions('fsolve','Display','off');
% options = odeset('RelTol',1e-3,'AbsTol',1e-6);

% Prepare to compute the non-dimensional potential drops
philim = 50;
Vrange = linspace(-philim,philim,1001);
Qrange = sign(Vrange).*sqrt(2).*sqrt(exp(Vrange)-1-Vrange);
V1fun = @(Qp) -interp1(Qrange,Vrange,-OmegaE*Qp,'spline');
V2fun = @(Qp) -interp1(Qrange,Vrange,-Qp,'spline');
V3fun = @(Qp) interp1(Qrange,Vrange,Qp,'spline');
V4fun = @(Qp) -interp1(Qrange,Vrange,-OmegaH*Qp,'spline');

% Define the derivative of the ionic charge density
dQdt = @(t,Q,psi) 2*psi(t)-V1fun(Q)-V2fun(Q)-V3fun(Q)-V4fun(Q);

if any(Vdc)
    % Find the steady-state charge density at the DC voltage
    Q = fsolve(@(Q) dQdt(0,Q,@(t) Vap2psi(Vdc)),0,optimopts);

else
    % Find the steady-state charge density at the initial voltage point
    Qinit = fsolve(@(Q) dQdt(0,Q,psi),0,optimopts);
    
    % Solve the charge evolution ODE to find Qp(t)
    sol = ode15s(@(t,Q) dQdt(t,Q,psi),time,Qinit);
    Q = deval(sol,time);
    
    if isnan(Q(end))
        fprintf('Not able to calculate Q(t) for the whole time range.\n')
        return
    end
end

% Calculate Debye layer potential drops from Q(t)
V1 = V1fun(Q);
V2 = V2fun(Q);
V3 = V3fun(Q);
V4 = V4fun(Q);

% Check that the potentials lie within the expected range
phimax = max(abs([V1,V2,V3,V4]));
if phimax > philim
    warning(['Need to increase philim above' num2str(phimax) '.']);
end

% Output the dimensional values
Q = sqrt(q*N0*epsp*VT)*Q; % charge density [C m-2]
Vs = struct('V1',V1*VT, 'V2',V2*VT, 'V3',V3*VT, 'V4',V4*VT);

end
