function Z = Z_from_J(sol)

nwaves = 2; % number of complete waves to analyse 

ind = (length(sol.J)-nwaves*100):length(sol.J);

J = sol.J(ind)*1e-3; % convert to units of Acm-2
V = sol.V(ind);
t = sol.time(ind)-sol.time(ind(1));

% --- get voltage wave parameters ---
% V = V0 + Vp*sin(omega*t)
V0 = sol.params.applied_voltage{3}; % voltage at the start of the first wave
Vp = sol.params.applied_voltage{end}-V0; % voltage amplitude
omega = 2*pi/sol.params.applied_voltage{end-1};

% --- get current wave parameters ---
% J = J0 + Jp*sin(omega*t + theta)
SineParams = sineFit(t,J');
J0 = SineParams(1);
Jp = SineParams(2);
omega = 2*pi*SineParams(3);
theta = SineParams(4)+pi; % + pi to account for the negative definition of J

Z = Vp/Jp*exp(-i*theta); % output in units of Ohm cm2

end
