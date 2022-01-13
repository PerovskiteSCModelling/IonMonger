function plot_IS(sol)
% A function to plot data from impedance spectroscopy simulations. The
% parameters of the Nyquist and Bode plots can easily be edited and new
% plots can be added below. `sol` is a solution structure from an IS
% simulation. This function works for full solutions or reduced solutions
% obtained using the `reduced_output` parameter.

% === extract data ===

% check which type of solution structure
if size(sol,2)>1
    % recieved full solution array
    [X,R] = impedance_analysis(sol);
    for j = 1:size(sol,2)
        freqs(j) = sol(j).freq;
    end
else
    if ~isfield(sol,'V')
        % recieved reduced solution structure
        [X,R,freqs] = struct2array(sol,{'X','R','freqs'});
    else
        % recieved non-impedance solution structure
        error('plot_IS was given a solution structure that was not from an impedance simulation.')
    end
end

% === default plots ===

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

T.TileSpacing = 'compact';

% Nyquist plot

figure(97)
plot(R,-X,'-sr')
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',18,'DataAspectRatio',[1 1 1])
ylabel('-X / $\Omega$cm$^2$','Interpreter','latex')
xlabel('R / $\Omega$cm$^2$','Interpreter','latex')

% === additional plots ===
% Users can create their own plots here. Any variables other than `X`, `R`,
% and `freqs` will need to be extracted from `sol`.


end