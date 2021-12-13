function fit = FourierFit(t,S)
% A function to fit a sinusoidal curve to a signal using the Fourier
% transform. Here `t` is a time vector and `S` is the signal vector. The
% function returns a structure, `fit` containing all the information about
% the fitted curve and the error. Note that this function is only intended
% to analyse signals where a signal frequency is dominant.

if size(S,1)>1
    % ensure S is a row vector
    S = S';
end
if size(t,1)>1
    % ensure t is a row vector
    t = t';
end

S0 = mean(S);
S = S-S0;

L = length(t);
T = t(2)-t(1); % sample period
Fs = 1/T; % sample frequency

Y = fft(S); % perform fast Fourier transform

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L; % frequency space

i = (P1/max(P1)>0.1); % indices of frequencies with 10% significance
[~,i] = maxk(P1,1);
Sp = P1(i); % amplitudes
omega = f(i); % frequencies
theta = angle(Y(i))+pi/2; % phase
theta(theta>pi) = theta(theta>pi)-2*pi;

X = 0;
for k = 1:length(omega)
    X = X + Sp(k)*sin(2*pi*omega(k)*t+theta(k));
end

% return information such that S = S0+Sp*sin(2*pi*omega*t+theta)
fit.S = @(t) S0+Sp*sin(2*pi*omega*t+theta);
fit.omega = omega; % frequency
fit.S0 = S0; % average signal value
fit.Sp = Sp; % signal amplitude
fit.theta = theta; % signal phase
S = S+S0;
fit.err = norm(S-fit.S(t),2)/norm(S); % fit error

end