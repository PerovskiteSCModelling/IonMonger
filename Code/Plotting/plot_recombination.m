function plot_recombination(sol)
% Plots the rates of different types of bulk recombination as functions of
% time and space.

% Check sol structure
if size(sol,2)>1 % received structure array from IS simulation
    error(['plot_recombination was given a solution structure array from an ' ...
        'impedance spectroscopy simulation. To use plot_recombination for the '...
        'n-th sample frequency solution, use `plot_recombination(sol(n),...)`'])
elseif isfield(sol,'X') % received reduced solution structure from IS simulation
    error(['plot_recombination was given a reduced solution structure from an ' ...
        'impedance spectroscopy simulation. To use plot_recombination with an ' ...
        'IS solution, ensure reduced_output=false'])
end

% Unpack the parameters and dimensional variables
[G0, n0, p0, N0, brate, ni2, gamma, tor, tor3, Cn, Cp, R, SRH, Auger, jay] = ...
	struct2array(sol.params, {'G0','n0','p0','N0','brate','ni2','gamma', ...
                              'tor','tor3','Cn','Cp','R','SRH','Auger','jay'});
x = sol.vectors.x;
[n, p, P] = struct2array(sol.dstrbns, {'n','p','P'});
time = sol.time;

% Convert to dimensionless variables (see the scaling in numericalsolver.m)
n = n/n0;
p = p/p0;
P = P/N0;

% Calculate dimensional rates of recombination (see nondimensionalise.m)
R_bim = G0*brate*(n.*p-ni2);
R_Aug = G0*Auger(n,p,Cn,Cp,ni2);
R_SRH = G0*SRH(n,p,gamma,ni2,tor,tor3);
R_tot = G0*R(n,p,P); % using dimensionless variables

% Rescale back to dimensional variables
% n = n0*n;
% p = p0*p;
% P = N0*P;

% Set default figure options
set(0,'defaultAxesFontSize',10); % Make axes labels smaller
set(0,'defaultTextInterpreter','latex'); % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex'); % For latex tick labels
set(0,'defaultLegendInterpreter','latex'); % For latex legends

% Plot the rates of recombination in space and time
figure;
subplot(2,2,1);
surf(x,time,R_tot,'LineStyle','none');
xlabel('Distance (nm)'); ylabel('Time (s)');
title('Total recombination rate');
subplot(2,2,2);
surf(x,time,R_bim,'LineStyle','none');
xlabel('Distance (nm)'); ylabel('Time (s)');
title('Bimolecular recombination rate');
subplot(2,2,3);
surf(x,time,R_Aug,'LineStyle','none');
xlabel('Distance (nm)'); ylabel('Time (s)');
title('Auger recombination rate');
subplot(2,2,4);
surf(x,time,R_SRH,'LineStyle','none');
xlabel('Distance (nm)'); ylabel('Time (s)');
title('SRH recombination rate');
sgtitle('Recombination rates (m$^{-3}\,$s$^{-1}$)','Interpreter','latex');

% Integrate to find the current lost to each type of recombination
J_tot = -jay*trapz(x/x(end),R_tot'/G0);
J_bim = -jay*trapz(x/x(end),R_bim'/G0);
J_Aug = -jay*trapz(x/x(end),R_Aug'/G0);
J_SRH = -jay*trapz(x/x(end),R_SRH'/G0);

% Reset default figure options
set(0,'defaultAxesFontSize',18); % Make axes labels larger

% Plot the integrated rates of recombination over time
figure;
hold on;
plot(time,J_tot,'DisplayName','Total');
plot(time,J_bim,'DisplayName','Bimolecular');
plot(time,J_Aug,'DisplayName','Auger');
plot(time,J_SRH,'DisplayName','SRH');
for i = 1:100:length(time) % section markers
    plot(time([i,i]),ylim,'k--','HandleVisibility','off');
end
ylabel('Current density (mA$\cdot$cm$^{-2}$)');
xlabel('Time (s)');
legend('Location','Best','FontSize',12);

end
