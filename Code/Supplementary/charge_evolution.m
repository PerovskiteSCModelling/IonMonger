function [Q, Vs] = charge_evolution(params)
% A function to compute the Debye layer charge density over time.

% Unpack parameters
[psi, time, OmegaE, OmegaH, VT] = ...
    struct2array(params,{'psi','time','OmegaE','OmegaH','VT'});

% Solver settings
optimopts = optimoptions('fsolve','Display','off');
% options = odeset('RelTol',1e-3,'AbsTol',1e-6);

% Find the steady-state charge density at the initial voltag point
Qinit = fsolve(@(Q) dQdt(0,Q,psi,OmegaE,OmegaH), 0, optimopts);

% Solve the charge evolution ODE to find Qp(t)
sol = ode15s(@(t,Q) dQdt(t,Q,psi,OmegaE,OmegaH),time,Qinit);
Q = deval(sol,time);

if isnan(Q(end))
    fprintf('Not able to calculate Q(t) for the whole time range.\n')
    return
end

% Calculate Debye layer potential drops from Q(t)
[V1, V2, V3, V4] = VofQ(Q,OmegaE,OmegaH);
Vs = struct('V1',V1*VT, 'V2',V2*VT, 'V3',V3*VT, 'V4',V4*VT);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions used by the code above:

function rhs = dQdt(t,Q,psi,OmegaE,OmegaH)
[V1, V2, V3, V4] = VofQ(Q,OmegaE,OmegaH);
rhs = 2*psi(t)-V1-V2-V3-V4;
end

function [V1, V2, V3, V4] = VofQ(Qp,OmegaE,OmegaH)
V = linspace(-50,50,101);
Q = sign(V).*sqrt(2).*sqrt(exp(V)-1-V);
V1 = -interp1(Q,V,-OmegaE*Qp,'spline');
V2 = -interp1(Q,V,-Qp,'spline');
V3 = interp1(Q,V,Qp,'spline');
V4 = -interp1(Q,V,-OmegaH*Qp,'spline');
end
