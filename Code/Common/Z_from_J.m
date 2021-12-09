function Z = Z_from_J(sol)
% This function extracts the impedance from the solution to a single
% impedance measurement. `sol` is the solution structure with a sinusoidal
% voltage protocol.

nwaves = 2; % number of complete waves to analyse 

ind = (length(sol.J)-nwaves*100):length(sol.J);     % indices of the timesteps
                                                    % to be analysed

J = sol.J(ind)*1e-3; % convert to units of Acm-2
V = sol.V(ind);
t = sol.time(ind)-sol.time(ind(1));

% --- get voltage wave parameters ---
% V = V0 + Vp*sin(omega*t)
V0 = sol.V(1); % voltage at the start of the first wave
Vp = sol.params.applied_voltage{end}-V0; % voltage amplitude
omega = 2*pi/sol.params.applied_voltage{end-1};

% --- get current wave parameters ---
% J = J0 + Jp*sin(omega*t + theta)

% Fourier transform
a0 = 1/t(end)*trapz(t,J');
a1 = 2/t(end)*trapz(t,J'.*cos(omega*t));
b1 = 2/t(end)*trapz(t,J'.*sin(omega*t));

J0 = a0;
Jp = sqrt(a1^2+b1^2);
if b1>0, theta = atan(a1/b1);
elseif a1>0, theta = atan(a1/b1)+pi;
elseif a1<0, theta = atan(a1/b1)-pi;
else, error('Current phase calculation was unsuccessful'), end

theta = theta + pi; % to account for negative definition of current

Z = Vp/Jp*exp(-i*theta); % output impedance in units of Ohm cm2

end