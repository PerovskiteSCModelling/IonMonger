function plot_nEIS(sol)
% A function to plot data from impedance spectroscopy simulations. `sol` is
% a solution structure from an IS simulation. This function works for full
% solutions or reduced solutions obtained using the `reduced_output`
% parameter.


%% Extract data

% Check which type of solution structure
if size(sol,2)>1
    % received full solution array
    [X,R] = impedance_analysis(sol);
    for j = 1:size(sol,2)
        freqs(j) = 1/sol(j).params.applied_voltage{2};
    end
else
    if ~isfield(sol,'V')
        % received reduced solution structure
        [X, R, freqs] = struct2array(sol,{'X','R','freqs'});
    else
        % received non-impedance solution structure
        error(['plot_nEIS was given a solution structure that was ', ...
               'not from an impedance simulation.'])
    end
end


%% Default plots

% Set default figure options
set(0,'defaultAxesFontSize',14); % Make axes labels larger
set(0,'defaultTextInterpreter','latex'); % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex'); % For latex tick labels
set(0,'defaultLegendInterpreter','latex'); % For latex legends
M = 2; % marker size
L = 0.5; % line width

% Set label names for the different orders
Rlabel = {'Re(Z$_{1}$) / $\Omega$cm$^2$', ...
          'Re(Z$_{2}$) / V$\Omega$cm$^2$', ...
          'Re(Z$_{3}$) / V$^2\Omega$cm$^2$'};
Xlabel = {'-Im(Z$_{1}$) / $\Omega$cm$^2$', ...
          '-Im(Z$_{2}$) / V$\Omega$cm$^2$', ...
          '-Im(Z$_{3}$) / V$^2\Omega$cm$^2$'};

% Choose line colour
colour = [126,47,142]/255; % purple
% colour = [0,0,0]/255; % black
style = 'o-';

for j = 1:size(X,2)
    figure('Name',['Order ' num2str(j)]);
    % Nyquist plot
    subplot(2,1,1); hold on;
    plot(R(:,j),-X(:,j),style,'LineWidth',L,'MarkerSize',M, ...
        'Color',colour,'MarkerFaceColor',colour);
    set(gca,'DataAspectRatio',[1 1 1]);
    ylabel(Xlabel{j});
    xlabel(Rlabel{j});
    % Frequency plots
    subplot(2,2,3); hold on;
    plot(freqs,R(:,j),style,'LineWidth',L,'MarkerSize',M, ...
        'Color',colour,'MarkerFaceColor',colour);
    xlabel('Frequency (Hz)');
    ylabel(Rlabel{j});
    set(gca,'Xscale','log');
    subplot(2,2,4); hold on;
    plot(freqs,-X(:,j),style,'LineWidth',L,'MarkerSize',M, ...
        'Color',colour,'MarkerFaceColor',colour);
    xlabel('Frequency (Hz)');
    ylabel(Xlabel{j});
    set(gca,'Xscale','log');
end

end
