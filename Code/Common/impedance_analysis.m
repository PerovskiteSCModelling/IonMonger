function [X,R] = impedance_analysis(sol,nwaves)
% This function calculates the impedance from multiple impedance
% measurements and decomposes it into its real and imaginary components.
% sol is a structure array containing the solutions to multiple impedance
% measurements where sol(i) is the solution structure of the i-th sample
% frequency. The function returns `X`, the imaginary component of impedance
% at each sample frequency, and `R`, the real component of impedance at
% each sample frequency, both in units of Ohm cm2. `nwaves` is an optional
% argument specifying the number of complete periods on which to perform
% the phase analysis. The default value is 2.

if nargin<2 ; nwaves = 2; end % number of complete periods to analyse

% Extract the impedance protocol
nf = sol(1).impedance_protocol{6}; % number of frequencies to be sampled
min_f = sol(1).impedance_protocol{2}; % minimum sample frequency
max_f = sol(1).impedance_protocol{3}; % maximum frequency
Vp = sol(1).impedance_protocol{5}; % voltage amplitude
V0 = sol(1).impedance_protocol{4}; % DC voltage

% Check number of periods
if nwaves > sol(1).impedance_protocol{7}
    error(['impedance_analysis was asked to analyse ' num2str(nwaves), ...
        ' complete periods but the solution only contains ', ...
        num2str(sol(1).impedance_protocol{7}) ' complete periods.'])
end

% Define the frequency range
freqs = logspace(log10(min_f),log10(max_f),nf);

% Compute the complex impedance
order = 3;
Z = nan(nf,order)+1i*nan(nf,order); % preallocate Z
for j = 1:nf
    try % Extract the impedance from each measurement
        
        % Get times and current (= -photocurrent)
        ind = (length(sol(j).J)-nwaves*100):length(sol(j).J);
        t = sol(j).time(ind)-sol(j).time(ind(1));
        J = -sol(j).J(ind)*1e-3; % [Acm-2]
        
        % Perform sinusoidal fit via Fourier transform
        fit = FourierFit(t,J,freqs(j));
        theta = fit.theta;
        Jp = fit.Sp; % extract sinusoidal current amplitude
        
        if fit.err>1e-1
            warning(['sinusoidal fit of current output for frequency '...
                num2str(j) ' may be inaccurate'])
        end
        
        % Output impedance
        Z(j,1) = Vp/Jp*exp(1i*theta); % [Ohm cm2]
        
        % Fit the second-order response
        fit2 = FourierFit(t,J-transpose(fit.S(t)),2*freqs(j));
        Jp2 = fit2.Sp;
        theta2 = fit2.theta;                       
        Z(j,2) = Vp/(Jp2*exp(-1i*theta2))^2;
        
        % Fit the third-order response
        fit3 = FourierFit(t,J-transpose(fit.S(t)+fit2.S(t)),3*freqs(j));
        Jp3 = fit3.Sp;
        theta3 = fit3.theta;
        Z(j,3) = Vp/(Jp3*exp(-1i*theta3))^3;
        
    catch me
        warning(['phase analysis of frequency ' num2str(j) ' was unsuccessful'])
            disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end
end

% Decompose into real and imaginary components
R = real(Z);
X = imag(Z);

end
