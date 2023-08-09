function [X,R] = impedance_analysis(sol,nwaves)
% This function calculates the nonlinear, complex impedances from multiple
% impedance measurements and decomposes them into their real and imaginary
% components. The input sol is a structure array containing the solutions
% to multiple impedance measurements where sol(i) is the solution structure
% of the i-th sample frequency. The outputs are: `X`, the imaginary
% component of the nonlinear impedances at each sample frequency, and `R`,
% the real component of nonlinear impedances at each sample frequency. The
% optional argument `nwaves` can be used to specify the number of periods
% on which to perform the phase analysis; the default value is 2.

if nargin<2 ; nwaves = 2; end % number of complete periods to analyse

% Extract values from the impedance protocol
nf = sol(1).impedance_protocol{6}; % number of frequencies to be sampled
nper = sol(1).impedance_protocol{7}; % number of periods

% Check number of periods
if nwaves > nper
    error(['Cannot analyse ' num2str(nwaves) ' periods because the ' ...
           'solution only contains ' num2str(nper) ' complete periods.'])
end

% Define the frequency range and preallocate Z
order = 3;
Z = nan(nf,order)+1i*nan(nf,order);

% Compute the complex impedance from each measurement
for j = 1:nf
    try
        % Define length of input
        L = nwaves*100;

        % Get times and current (= negative photocurrent)
        ind = (length(sol(j).J)-L):(length(sol(j).J)-1);
        J = -sol(j).J(ind)*1e-3; % [Acm-2]
        V = sol(j).V(ind);       % [V]
        
        % Compute the Fourier coefficients 
        Jfft = fft(J)/L;
        Vfft = fft(V)/L;
        
        % Output the first, second and third-order impedance
        Z(j,1) = Vfft(nwaves+1)/(Jfft(nwaves+1));      % [Ohm cm2]
        Z(j,2) = Vfft(nwaves+1)^2./(Jfft(2*nwaves+1)); % [Ohm V cm2]
        Z(j,3) = Vfft(nwaves+1)^3./(Jfft(3*nwaves+1)); % [Ohm V2 cm2]
        
    catch me
        warning(['Phase analysis of frequency ' num2str(j) ' was unsuccessful.']);
        disp(getReport(me, 'extended', 'hyperlinks', 'on'));
    end
end

% Decompose into real and imaginary components
R = real(Z);
X = imag(Z);

end
