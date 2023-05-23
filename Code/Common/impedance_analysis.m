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
        ind = (length(sol(j).J)-nwaves*100):length(sol(j).J);
        t = sol(j).time(ind)-sol(j).time(ind(1));
        J = -sol(j).J(ind)*1e-3; % [Acm-2]
        V = sol(j).V(ind);

        % Use fft to get frequency information
        L = nwaves*100; 
        Jfft = 2*fft(J(1:end-1))/L;
        Vfft = 2*fft(V(1:end-1))/L;
        
        
        % Output impedance
        Z(j,1) = Vfft(3)/(Jfft(3)); % [Ohm cm2]
        
                  
        Z(j,2) = Vfft(3)^2./(Jfft(5));
        

        Z(j,3) = Vfft(3)^3./(Jfft(7));
        
%         if j==1.5*64
%             % Plot harmonic response using the discrete Fourier transform
%             FTJ = fft(J-mean(J));
%             P2 = abs(FTJ/length(J)); % two-sided spectrum
%             P1 = P2(1:floor(length(J)/2)+1);
%             P1(2:end-1) = 2*P1(2:end-1); % one-sided spectrum
%             P1 = P1*1e3; % rescale to mAcm-2 for plotting
%             Fs = 1/t(2); % sampling frequency
%             f = Fs*(0:floor(length(J)/2))/length(J);
%             figure; hold on;
%             cut = 20;
%             plot(f(1:cut),P1(1:cut));
%             plot(freqs(j)*[1,1],[0,max(P1)],'--'); % fundamental frequency
% %             title('Single-sided amplitude spectrum of J(t)');
%             ylabel('Amplitude (mAcm$^{-2}$)');
%             xlabel('Frequency (Hz)');
%             legend('Amplitude spectrum of J(t)','Fundamental frequency');
%             drawnow;
%         end
        
    catch me
        warning(['phase analysis of frequency ' num2str(j) ' was unsuccessful'])
            disp( getReport( me, 'extended', 'hyperlinks', 'on' ) )
    end
end

% Decompose into real and imaginary components
R = real(Z);
X = imag(Z);

end
