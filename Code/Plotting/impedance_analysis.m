function [X,R] = impedance_analysis(sol)
% This function calculates the impedance from multiple impedance
% measurements, decomposes it into its real and imaginary components and
% plots the results in the form of Nyquist and Bode plots. `sol` is a
% structure array containing the solutions to multiple impedance
% measurements where `sol(i)` is the solution structure of the i-th sample
% frequency. The functionreturns `X`, the imaginary component of impedance
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
        nwaves = 2; % number of complete waves to analyse 
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
        % do nothing
    end
end

% decompose into real and imaginary components
R = real(Z);
X = imag(Z);

% === Plot the data ===

% Bode plot
figure(96)
T = tiledlayout(2,1);
ax1 = nexttile;
plot(ax1,freqs,X,'x-b')
ax2 = nexttile;
plot(freqs,R,'x-b')

set(ax1,'XScale','log','TickLabelInterpreter','latex','FontSize',18,'YDir','reverse')
set(ax2,'XScale','log','TickLabelInterpreter','latex','FontSize',18)

ylabel(ax1,'X / $\Omega$cm$^2$','Interpreter','latex')
ylabel(ax2,'R / $\Omega$cm$^2$','Interpreter','latex')
xlabel(ax2,'frequency / Hz','Interpreter','latex')
ylim(ax2,[0 inf])
title(ax1,['$V_{DC} =$ ' num2str(V0) 'V'],'Interpreter','latex')

T.TileSpacing = 'compact';
drawnow;

% Nyquist plot

figure(97)
plot(R,-X,'-sr')
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',18)
ylabel('-X / $\Omega$cm$^2$','Interpreter','latex')
xlabel('R / $\Omega$cm$^2$','Interpreter','latex')
title(['$V_{DC} =$ ' num2str(V0) 'V'],'Interpreter','latex')
drawnow;

end