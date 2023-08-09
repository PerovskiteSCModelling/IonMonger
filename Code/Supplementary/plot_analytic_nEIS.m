function plot_analytic_nEIS(params,figure_handles)
% A function to plot the first and second-order impedance spectra using
% analytic expressions.


%% Define the frequency range

% Define the sample frequencies, logarithmically spaced
freqs = logspace(log10(1e-4),log10(1e7),100); % [s-1]

% Define the corresponding angular frequencies
omega = 2*pi*freqs;


%% Compute the characteristic parameters

% Select the dominant type of recombination (Rb, Rp, Rn, Rl, or Rr)
rtype = 'Rl';

% Set the recombination parameters
[Jd, Fi, nid] = recombination_type(rtype,params);

% Compute the characteristic parameters
% Unpack parameters
[time, psi, psi2Vap, VT, b, q, epsp, DI, N0] = ...
    struct2array(params, {'time','psi','psi2Vap','VT','b','q','epsp', ...
                          'DI','N0'});

% Define the steady state voltage and half the perturbation amplitude
Vdc = psi2Vap(psi(time(end)));
Vp = max(psi2Vap(psi(time(end-100:end)))-Vdc)/2;

% Compute the steady-state ionic charge density and Debye layer potentials
params.Vdc = Vdc;
[Qdc, Vs] = ionic_charge(params);

% Compute the first and second derivatives of the potentials
params.Vdc = Vdc+Vp;
[Qdc2, Vs2] = ionic_charge(params);
params.Vdc = Vdc-Vp;
[Qdc1, Vs1] = ionic_charge(params);
FT = @(Vs) Vs.V1+Vs.V2+Vs.V3+Vs.V4;
deriv = @(F,V2,V1,Q2,Q1) (F(V2)-F(V1))./(Q2-Q1);
FTd = (deriv(FT,Vs2,Vs,Qdc2,Qdc)+deriv(FT,Vs,Vs1,Qdc,Qdc1))/2;
Fid = (deriv(Fi,Vs2,Vs,Qdc2,Qdc)+deriv(Fi,Vs,Vs1,Qdc,Qdc1))/2;
FTdd = (deriv(FT,Vs2,Vs,Qdc2,Qdc)-deriv(FT,Vs,Vs1,Qdc,Qdc1))./((Qdc2-Qdc1)/2);
Fidd = (deriv(Fi,Vs2,Vs,Qdc2,Qdc)-deriv(Fi,Vs,Vs1,Qdc,Qdc1))./((Qdc2-Qdc1)/2);

% Compute the model parameters
% nec = FT(Vs)./Fi(Vs);           % ectypal factor [non-dim.]
nap = FTd./Fid;                 % measured ectypal factor [non-dim.]
Jr0 = Jd.*exp(-Fi(Vs)/VT)*1e-3; % dark current density [A cm-2]
Gd = q*DI*N0/(VT*b);            % ionic conductance per unit area [A V-1 m-2]

% Compute the linear characteristic parameters
RLf = VT./Jr0.*(nap-nid);       % low-frequency resistance [Ohm cm2]
RHf = nid*VT./Jr0;              % high-frequency resistance [Ohm cm2]
TLf = 1./(Gd*nid*Fid);          % low-frequency time constant [s]
THf = epsp*VT*nid./(b*Jr0*1e4); % high-frequency time constant [s]
RT = RLf+RHf;                   % total resistance [Ohm cm2]

% Compute the nonlinear characteristic parameters
S1 = 1-Fidd.*(Gd*Jr0.*RHf.*TLf).^2/VT; % [non-dim.]
S2 = FTdd.*(Gd*RHf.*TLf).^2.*Jr0;      % [Ohm cm2]


%% Compute the impedance

% Compute the linear impedance
beta = TLf*omega;
gamma = 1-THf*TLf*omega.^2;
R(:,1) = (RHf*beta.^2+RT*gamma)./(beta.^2+gamma.^2); % [Ohm cm2]
X(:,1) = -beta.*(RT-RHf*gamma)./(beta.^2+gamma.^2);  % [Ohm cm2]

% Compute the second-order impedance
D1 = RT*S1+S2-(RT+4*RHf+4*S2*THf/TLf)*(TLf*omega).^2;
% D1 = RT*S1+S2-(RT+4*RHf)*beta.^2;
D2 = RT+RHf*S1+S2-RHf*beta.^2;
Z2 = 2*Jr0*(RT+1i*RHf*beta).^2.*(RT+2i*RHf*beta) ...
        .*(D1-2i*D2.*beta)./(D1.^2+(2*D2.*beta).^2);
R(:,2) = real(Z2);
X(:,2) = imag(Z2);


%% Default plots

% Set default figure options
set(0,'defaultAxesFontSize',14); % Make axes labels larger
set(0,'defaultTextInterpreter','latex'); % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex'); % For latex tick labels
set(0,'defaultLegendInterpreter','latex'); % For latex legends
M = 2; % marker size
L = 0.5; % line width

if nargin<2
    % Set up figure handles if they do not exist already
    Rlabel = {'Re(Z$_{1}$) / $\Omega$cm$^2$', ...
              'Re(Z$_{2}$) / V$\Omega$cm$^2$', ...
              'Re(Z$_{3}$) / V$^2\Omega$cm$^2$'};
    Xlabel = {'-Im(Z$_{1}$) / $\Omega$cm$^2$', ...
              '-Im(Z$_{2}$) / V$\Omega$cm$^2$', ...
              '-Im(Z$_{3}$) / V$^2\Omega$cm$^2$'};
    
    for j = 1:2
        figure('Name',['Order ' num2str(j)]);
        figure_handles(j) = gcf().Number;
        % Nyquist plot
        subplot(2,1,1); hold on;
        ylabel(Xlabel{j});
        xlabel(Rlabel{j});
        % Frequency plots
        subplot(2,2,3); hold on;
        xlabel('Frequency (Hz)');
        ylabel(Rlabel{j});
        subplot(2,2,4); hold on;
        xlabel('Frequency (Hz)');
        ylabel(Xlabel{j});
    end
end

% Plot the impedances
for j = 1:2
    % Nyquist plot
    figure(figure_handles(j));
    subplot(2,1,1); hold on;
    plot(R(:,j),-X(:,j),'-xb','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','b');
    set(gca,'DataAspectRatio',[1 1 1]);
    subplot(2,2,3); hold on;
    plot(freqs,R(:,j),'-xb','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','b');
    set(gca,'Xscale','log');
    subplot(2,2,4); hold on;
    plot(freqs,-X(:,j),'-xb','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','b');
    set(gca,'Xscale','log');
end

end

