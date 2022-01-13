function [X,R] = impedance_analysis(sol)
% This function calculates the impedance from multiple impedance
% measurements and decomposes it into its real and imaginary components. `sol` is a
% structure array containing the solutions to multiple impedance
% measurements where `sol(i)` is the solution structure of the i-th sample
% frequency. The function returns `X`, the imaginary component of impedance
% at each sample frequency, and `R`, the real component of impedance at
% each sample frequency.

% extract the impedance protocol
nf = sol(1).impedance_protocol{6}; % number of frequencies to be sampled
min_f = sol(1).impedance_protocol{2}; % minimum sample frequency
max_f = sol(1).impedance_protocol{3}; % maximum frequency
Vp = sol(1).impedance_protocol{5}; % voltage amplitude
V0 = sol(1).impedance_protocol{4}; % DC voltage
n_wave = sol(1).impedance_protocol{7}; % number of complete sine waves

freqs = logspace(log10(min_f),log10(max_f),nf);

Z = nan(nf,1)+i*nan(nf,1); % preallocate Z
for j = 1:nf
    try % extract the impedance from each measurement
        nwaves = 2; % number of complete periods to analyse 
        ind = (length(sol(j).J)-nwaves*100):length(sol(j).J); % indices of the timesteps
                                                              % to be analysed
        J = sol(j).J(ind)*1e-3; % convert to units of Acm-2
        t = sol(j).time(ind)-sol(j).time(ind(1));
        
        fit = FourierFit(t,J,freqs(j)); % perform sinusoidal fit via Fourier transform
        theta = fit.theta+pi; % add pi to account for negative current definition
        Jp = fit.Sp; % extract sinusoidal current amplitude
        
        if fit.err>1e-1
            warning(['sinusoidal fit of current output for frequency '...
                num2str(j) ' may be inaccurate'])
        end
        
        Z(j) = Vp/Jp*exp(-i*theta); % output impedance in units of Ohm cm2
    catch me
        warning(['phase analysis of frequency ' num2str(j) ' was unsuccessful'])
            disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end
end

% decompose into real and imaginary components
R = real(Z);
X = imag(Z);

end