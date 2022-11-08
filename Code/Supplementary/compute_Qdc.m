function [Qdc, Vs] = compute_Qdc(Vap,params)
% A function to compute the steady-state ionic charge density in the Debye
% layers, as well as the potential drops V1-4, by solving the asymptotically
% reduced problem and neglecting series resistance.

% Parameter input
[Vap2psi, OmegaE, OmegaH, q, LD, N0, VT] = ...
    struct2array(params,{'Vap2psi','OmegaE','OmegaH','q','LD','N0','VT'});

% Compute the vector of non-dimensional voltage points
psit = Vap2psi(Vap);

% Program settings
options = optimoptions(@fsolve,'Display','off');

% Prepare to compute the non-dimensional potential drops
Vrange = linspace(-50,50,101);
Qrange = sign(Vrange).*sqrt(2).*sqrt(exp(Vrange)-1-Vrange);
Vpfun = @(Qp) interp1(Qrange,Vrange,Qp,'spline');
Vmfun = @(Qp) interp1(Qrange,Vrange,-Qp,'spline');
VEfun = @(Qp) -interp1(Qrange,Vrange,-OmegaE*Qp,'spline');
VHfun = @(Qp) interp1(Qrange,Vrange,-OmegaH*Qp,'spline');

% Preallocation vectors
[Qdc, Vm, Vp, VE, VH] = deal(NaN(size(psit)));

% Solve for the steady-state ionic charge density at each voltage point
Qinit = 0;
for i = 1:length(psit)
    Qdc(i) = fsolve(@(Q) 2*psit(i)+Vmfun(Q)-Vpfun(Q)-VEfun(Q)+VHfun(Q), ...
                    Qinit, options);
    Qinit = Qdc(i);
    Vm(i) = Vmfun(Qinit);
    Vp(i) = Vpfun(Qinit);
    VE(i) = VEfun(Qinit);
    VH(i) = VHfun(Qinit);
end

% Re-dimensionalise
Qdc = -q*LD*N0*1e3*Qdc; % check scaling
Vs = struct('V1',VT*VE, 'V2',-VT*Vm, 'V3',VT*Vp, 'V4',-VT*VH);

end
