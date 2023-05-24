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

% Define the frequency range and preallocate Z
freqs = logspace(log10(min_f),log10(max_f),nf);
order = 3;
Z = nan(nf,order)+1i*nan(nf,order);

% Compute the complex impedance
for j = 1:nf
    try % Extract the impedance from each measurement
        
        % Get times and current (= -photocurrent)
        ind = (length(sol(j).J)-nwaves*100):(length(sol(j).J)-1);
        J = -sol(j).J(ind)*1e-3; % [Acm-2]
        V = sol(j).V(ind);

        % Use fft to get frequency information
        L = nwaves*100; 
        Jfft = 2*fft(J)/L;
        Vfft = 2*fft(V)/L;
        
        % Output first, second and third-order impedance
        Z(j,1) = Vfft(nwaves+1)/(Jfft(nwaves+1)); % [Ohm cm2]     
        Z(j,2) = Vfft(nwaves+1)^2./(Jfft(2*nwaves+1));
        Z(j,3) = Vfft(nwaves+1)^3./(Jfft(3*nwaves+1));
        
    catch me
        warning(['phase analysis of frequency ' num2str(j) ' was unsuccessful'])
            disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end
end

% Decompose into real and imaginary components
R = real(Z);
X = imag(Z);

end
