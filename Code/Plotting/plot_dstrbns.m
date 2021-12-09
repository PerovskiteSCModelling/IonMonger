function plot_dstrbns(sol,plotindex)
% Plots the solution variables as functions of space at the times specified
% by plotindex, which must be a vector of integers corresponding to the
% desired points in sol.time (for example, plotindex = [1,101,201,301]),
% and the current density and interfacial recombination losses over time.

% Retrieve the nondimensional parameter values
[kE, kH] = struct2array(sol.params, {'kE','kH'});

% Unpack the spatial vectors and dimensional solution variables
[x, xE, xH] = struct2array(sol.vectors,{'x','xE','xH'});
[P, phi, n, p, phiE, nE, phiH, pH] ...
    = struct2array(sol.dstrbns,{'P','phi','n','p','phiE','nE','phiH','pH'});
[time, J, Jl, Jr] = struct2array(sol,{'time','J','Jl','Jr'});

% Set default figure options
set(0,'defaultAxesFontSize',10); % Make axes labels smaller
set(0,'defaultTextInterpreter','latex') % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex') % For latex tick labels
set(0,'defaultLegendInterpreter','latex') % For latex legends

% Shading
if length(plotindex)>1
    shade = @(tt) double((tt-plotindex(1))/(plotindex(end)-plotindex(1)));
else
    shade = @(tt) 1;
end
P_colour = [1 0 1];
phi_colour = [119 172 48]/255; % green
n_colour = [0 0 1];
p_colour = [1 0 0];

% Plot the solution variables in space at the chosen times
figure;
subplot(2,2,1); hold on;
for tt = plotindex
    semilogy(x,P(tt,:),'color',shade(tt)*P_colour);
    set(gca,'YScale','log')
end
xlim([xE(1),xH(end)]);
xlabel('Distance (nm)'); ylabel('Ion Vacancy Density (cm$^{-3}$)');
subplot(2,2,2); hold on;
for tt = plotindex
    plot(x,phi(tt,:),'color',shade(tt)*phi_colour);
    plot(xE,phiE(tt,:),'color',shade(tt)*phi_colour);
    plot([x(end); xH],[phi(tt,end) phiH(tt,:)],'color',shade(tt)*phi_colour);
end
xlim([xE(1),xH(end)]);
xlabel('Distance (nm)'); ylabel('Electric Potential (V)');
subplot(2,2,3); hold on;
for tt = plotindex
    semilogy(x,n(tt,:),'color',shade(tt)*n_colour);
    semilogy(xE,nE(tt,:),'color',shade(tt)*n_colour);
    set(gca,'YScale','log')
end
xlim([xE(1),xH(end)]);
xlabel('Distance (nm)'); ylabel('Electron Concentration (cm$^{-3}$)');
subplot(2,2,4); hold on;
for tt = plotindex
    semilogy(x,p(tt,:),'color',shade(tt)*p_colour);
    semilogy(xH,pH(tt,:),'color',shade(tt)*p_colour);
    set(gca,'YScale','log')
end
xlim([xE(1),xH(end)]);
xlabel('Distance (nm)'); ylabel('Hole Concentration (cm$^{-3}$)');

% Reset default figure options
set(0,'defaultAxesFontSize',18); % Make axes labels larger

% Plot the current density and interfacial recombination losses over time
figure;
hold on;
plot(time,J,'k');
plot(time,Jl,'b--');
plot(time,Jr,'r--');
xlabel('Time (s)'); ylabel('Current Density (mAcm$^{-2}$)');
legend('J','-Rl','-Rr');

% Add markers to the current showing the time points chosen for plotting
plot(time(plotindex),J(plotindex),'mx','HandleVisibility','off');
for i = find(plotindex==1) % make up for missing J(1) with approximation
    plot(time(1),J(2),'mx','HandleVisibility','off');
end

end
