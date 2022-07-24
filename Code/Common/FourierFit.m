function fit = FourierFit(t,S,f)
% A function to fit a sinusoidal curve to a signal using the Fourier
% transform. Here `t` is a time vector and `S` is the signal vector. The
% function returns a structure, `fit` containing all the information about
% the fitted curve and the error. Note that this function is only intended
% to analyse signals where a signal frequency, `f`, is present. The 
% signal must contain an integer number of complete periods and the time 
% vector t must start from time = 0. In IonMonger's output, each period
% comprises 100 time steps, meaning N complete waves should contain 
% 100*N+1 points.

% This function obtains a fit of the form S = S0+Sp*sin(2*pi*f*t+theta).

% ===================== example of periodic analysis ======================
% To analyse the 'impedance' of any periodic variable in the solution
% structure, use and adapt the following code. Ensure `reduced_output` is
% set to `false` in the parameters file.

% % Set default figure options
% set(0,'defaultAxesFontSize',18); % Make axes labels larger
% set(0,'defaultTextInterpreter','latex'); % For latex axis labels
% set(0,'defaultAxesTickLabelInterpreter','latex'); % For latex tick labels
% set(0,'defaultLegendInterpreter','latex'); % For latex legends
% 
% nwaves = 2; % number of waves to analyse
% [f,Sp,S0,theta_S] = deal(NaN(size(sol))); % preallocate vectors
% 
% for j = 1:length(sol)
%     % set the variable you wish to extract phase and amplitude from
%     var = sol(j).J;     % e.g sol(j).dstrbns.P(:,end); or sol(j).Jd;
%                            
%     % gather indices of vector corresponding to periods to be analysed
%     ind = (length(var)-nwaves*100):length(var);
%     S = var(ind);    % construct signal S(t) with n waves
%      
%     % construct corresponding time vector
%     t = sol(j).time(ind)-sol(j).time(ind(1));  % ensures time starts at 0
%     
%     f(j) = 1/sol(j).params.applied_voltage{end-1}; % get frequency of input
%     fit = FourierFit(t,S,f(j)); % fit the signal
%     Sp(j) = fit.Sp;  
%     S0(j) = fit.S0;
%     theta_S(j) = fit.theta+pi; % get phase from fit
%     % (some quantities may require an extra phase offset of pi)
%     Vp = sol(1).impedance_protocol{5}; % get voltage amplitude
%     Z(j) = Vp/fit.Sp*exp(-1i*theta_S(j)); % create 'impedance'
% end
% 
% figure();    % Nyquist plot
% plot(real(Z),-imag(Z),'-sr');
% grid on;
% set(gca,'DataAspectRatio',[1 1 1]); % comment if poor view
% ylabel('-imag($V/S$)');
% xlabel('real($V/S$)');
% 
% figure();    % Real and imaginary frequency plot
% T = tiledlayout(2,1);
% ax1 = nexttile;
% semilogx(ax1,f,-imag(Z),'x-b');
% ax2 = nexttile;
% semilogx(ax2,f,real(Z),'x-b');
% ylabel(ax1,'-imag($V/S$)');
% ylabel(ax2,'real($V/S$)');
% xlabel(ax2,'frequency / Hz');
% 
% figure();    % Amplitude and phase frequency plot
% T = tiledlayout(2,1);
% ax1 = nexttile;
% semilogx(ax1,f,Sp,'x-b');
% ax2 = nexttile;
% semilogx(ax2,f,theta_S,'x-b');
% ylabel(ax1,'amplitude $S_p$');
% ylabel(ax2,'phase $\theta_S$');
% xlabel(ax2,'frequency / Hz');

% =========================================================================

if size(S,1)>1
    S = S'; % make S a row vector
end
if size(t,1)>1
    t = t'; % make t a row vector
end

a0 = 1/t(end)*trapz(t,S);
a1 = 2/t(end)*trapz(t,S.*cos(2*pi*f*t));
b1 = 2/t(end)*trapz(t,S.*sin(2*pi*f*t));

S0 = a0;
Sp = sqrt(a1^2+b1^2);
if b1>0, theta = atan(a1/b1);
elseif a1>0, theta = atan(a1/b1)+pi;
elseif a1<0, theta = atan(a1/b1)-pi;
else, error('Phase calculation was unsuccessful');
end

fit.S = @(t) S0+Sp*sin(2*pi*f*t+theta);
fit.omega = f; % frequency
fit.S0 = S0; % average signal value
fit.Sp = Sp; % signal amplitude
fit.theta = theta; % signal phase
fit.err = norm(abs(S-fit.S(t)),2)/norm(abs(S)); % fit error

end