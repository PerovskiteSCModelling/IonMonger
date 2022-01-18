function fit = FourierFit(t,S,omega)
% A function to fit a sinusoidal curve to a signal using the Fourier
% transform. Here `t` is a time vector and `S` is the signal vector. The
% function returns a structure, `fit` containing all the information about
% the fitted curve and the error. Note that this function is only intended
% to analyse signals where a signal frequency, `omega`, is present. The
% signal must contain an integer number of complete periods. In IonMonger's
% output, each period comprises 100 time steps, meaning N complete waves
% should contain 100*N+1 points.

% obtains a fit of the form S = S0+Sp*sin(2*pi*omega*t+theta)

% ===================== example of periodic analysis ======================
% To analyse the 'impedance' of any periodic variable in the solution
% structure, use and adapt the following code. Ensure `reduced_output` is
% set to `false` in the parameters file.
% 
% n = 2; % number of waves to analyse
% 
% for j = 1:length(sol)
%     % isolate a single variable as a function of time
%     var = sol(j).dstrbns.phi(:,end);
%     S = var(end-n*100:end); % construct signal from n waves
%     t = sol(j).time(end-n*100:end); % construct corresponding time vector
%     omega(j) = sol(j).freq; % get frequency of input
%     fit = FourierFit(t,S,omega(j)); % fit the signal
%     theta(j) = fit.theta; % get phase from fit (some quantities may require an extra phase offset of pi)
%     
%     Vp = sol(1).impedance_protocol{5}; % get voltage amplitude
%     Z(j) = Vp/fit.Sp*exp(-i*theta(j)); % create 'impedance'
% end
% 
% % Nyquist plot
% 
% figure()
% plot(real(Z),-imag(Z),'-sr')
% grid on
% ylabel('-X')
% xlabel('R')
% 
% % Bode plot
% figure()
% T = tiledlayout(2,1);
% ax1 = nexttile;
% semilogx(ax1,omega,-imag(Z),'x-b')
% ax2 = nexttile;
% semilogx(ax2,omega,real(Z),'x-b')
% ylabel(ax1,'-X')
% ylabel(ax2,'R')
% xlabel(ax2,'frequency / Hz')
% 
% =========================================================================

if size(S,1)>1
    % ensure S is a row vector
    S = S';
end
if size(t,1)>1
    % ensure t is a row vector
    t = t';
end

a0 = 1/t(end)*trapz(t,S);
a1 = 2/t(end)*trapz(t,S.*cos(2*pi*omega*t));
b1 = 2/t(end)*trapz(t,S.*sin(2*pi*omega*t));

S0 = a0;
Sp = sqrt(a1^2+b1^2);
if b1>0, theta = atan(a1/b1);
elseif a1>0, theta = atan(a1/b1)+pi;
elseif a1<0, theta = atan(a1/b1)-pi;
else, error('Phase calculation was unsuccessful'), end

fit.S = @(t) S0+Sp*sin(2*pi*omega*t+theta);
fit.omega = omega; % frequency
fit.S0 = S0; % average signal value
fit.Sp = Sp; % signal amplitude
fit.theta = theta; % signal phase
fit.err = norm(abs(S-fit.S(t)),2)/norm(abs(S)); % fit error

end