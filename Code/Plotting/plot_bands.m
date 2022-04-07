function plot_bands(sol,plotindex)
% Plots the energy levels as functions of space at the times specified
% by plotindex, which must be a vector of integers corresponding to the
% desired points in sol.time (for example, plotindex = [1,101,201,301]).

% Unpack the parameters, spatial vectors and dimensional solution variables
[EcE, Ec, Ev, EvH] = struct2array(sol.params,{'EcE','Ec','Ev','EvH'});
[x, xE, xH] = struct2array(sol.vectors,{'x','xE','xH'});
[phi, phiE, phiH, Efn, Efp, EfnE, EfpH] = ...
    struct2array(sol.dstrbns,{'phi','phiE','phiH','Efn','Efp','EfnE','EfpH'});

if ~any(Efn)
    % Compute and then load the quasi-Fermi levels
    sol = compute_QFLs(sol);
    [Efn, Efp, EfnE, EfpH] = ...
        struct2array(sol.dstrbns,{'Efn','Efp','EfnE','EfpH'});
end

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

% Plot the energy levels in space at the chosen times
figure;
hold on;
for tt = plotindex
    % Choose the reference energy level
    shift = phiE(tt,1); % sets the metal/ETL contact vacuum level as zero
    % Compute the vacuum level
    Evac  = shift-phi(tt,:);
    EvacE = shift-phiE(tt,:);
    EvacH = shift-phiH(tt,:);
    % Plot the vacuum level
    plot(x, Evac, '--','color',shade(tt)*phi_colour);
    plot(xE,EvacE,'--','color',shade(tt)*phi_colour);
    plot(xH,EvacH,'--','color',shade(tt)*phi_colour);
    % Plot the band energy levels
    plot(x, Evac+Ec,  'color',shade(tt)*phi_colour);
    plot(x, Evac+Ev,  'color',shade(tt)*phi_colour);
    plot(xE,EvacE+EcE,'color',shade(tt)*phi_colour);
    plot(xH,EvacH+EvH,'color',shade(tt)*phi_colour);
    % Plot the quasi-Fermi levels
    plot(x, shift+Efn(tt,:), 'color',shade(tt)*n_colour);
    plot(x, shift+Efp(tt,:), 'color',shade(tt)*p_colour);
    plot(xE,shift+EfnE(tt,:),'color',shade(tt)*n_colour);
    plot(xH,shift+EfpH(tt,:),'color',shade(tt)*p_colour);
end
xlabel('Distance (nm)'); ylabel('Energy (eV)');

% Reset default figure options
set(0,'defaultAxesFontSize',18); % Make axes labels larger

end
