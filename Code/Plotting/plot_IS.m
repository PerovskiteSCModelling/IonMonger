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
        freqs(j) = 1/sol(j).params.applied_voltage{2};
    end
else
    if ~isfield(sol,'V')
        % recieved reduced solution structure
        [X,R,freqs] = struct2array(sol,{'X','R','freqs'});
    else
        % recieved non-impedance solution structure
        error(['plot_IS was given a solution structure that was not from ',...
            'an impedance simulation.'])
    end
end

% === default plots ===
M = 2; % marker size
L = 0.5; % line width

% Bode plot
figure('Name','Bode plot')
T = tiledlayout(2,1);
ax1 = nexttile;
col = [0.3,0.7,0.1];
plot(ax1,freqs,abs(R+i*X),'-o','LineWidth',L,'MarkerSize',M,'MarkerFaceColor',col,'Color',col)
ax2 = nexttile;
plot(freqs,angle(R+i*X)/pi,'-o','LineWidth',L,'MarkerSize',M,'MarkerFaceColor',col,'Color',col)
set(ax1,'XScale','log','TickLabelInterpreter','latex','FontSize',18)
set(ax2,'XScale','log','TickLabelInterpreter','latex','FontSize',18)
ylabel(ax1,'$|Z|$ / $\Omega$cm$^2$','Interpreter','latex')
ylabel(ax2,'arg$(Z)$ / $\pi$','Interpreter','latex')
xlabel(ax2,'frequency / Hz','Interpreter','latex')
T.TileSpacing = 'compact';

% frequency plot
figure('Name','Frequency plot')
T = tiledlayout(2,1);
ax1 = nexttile;
plot(ax1,freqs,X,'-ob','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','b')
ax2 = nexttile;
plot(freqs,R,'-ob','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','b')
set(ax1,'XScale','log','TickLabelInterpreter','latex','FontSize',18,...
    'YDir','reverse')
set(ax2,'XScale','log','TickLabelInterpreter','latex','FontSize',18)
ylabel(ax1,'X / $\Omega$cm$^2$','Interpreter','latex')
ylabel(ax2,'R / $\Omega$cm$^2$','Interpreter','latex')
xlabel(ax2,'frequency / Hz','Interpreter','latex')
ylim(ax2,[0 inf])
T.TileSpacing = 'compact';

% Nyquist plot
figure('Name','Nyquist plot')
plot(R,X,'-or','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','r')
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',18,'DataAspectRatio',[1 1 1],'YDir','reverse')
ylabel('X / $\Omega$cm$^2$','Interpreter','latex')
xlabel('R / $\Omega$cm$^2$','Interpreter','latex')

% 3D plot
figure('Name','3D impedance plot')
plot3(freqs,R,X,'-o','LineWidth',L,'MarkerSize',M,'MarkerFaceColor','m','Color','m')
hold on
ax = gca;
plot3(ax.XLim(2)*ones(size(R)),R,X,'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto R-X plane
plot3(freqs,ax.YLim(1)*ones(size(R)),X,'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto f-X plane
plot3(freqs,R,ax.ZLim(2)*ones(size(R)),'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto f-R plane
patch([min(freqs) max(freqs) max(freqs) min(freqs)],min(R)*[1,1,1,1],[ax.ZLim(2),ax.ZLim(2),ax.ZLim(1),ax.ZLim(1)],...
    'r','FaceAlpha',0.1,'EdgeColor','none')
patch(max(freqs)*[1,1,1,1],[ax.YLim(1),ax.YLim(2),ax.YLim(2),ax.YLim(1)],[ax.ZLim(2),ax.ZLim(2),ax.ZLim(1),ax.ZLim(1)],...
    'b','FaceAlpha',0.1,'EdgeColor','none')
patch([min(freqs) max(freqs) max(freqs) min(freqs)],[ax.YLim(1),ax.YLim(1),ax.YLim(2),ax.YLim(2)],max(X)*[1,1,1,1],...
    'g','FaceAlpha',0.1,'EdgeColor','none')
set(gcf,'Renderer','painters')
hold off
set(gca,'XScale','log','ZDir','reverse','YDir','reverse','TickLabelInterpreter','latex','FontSize',13)
grid on
box on
zlabel('X / $\Omega$cm$^2$','Interpreter','latex')
ylabel('R / $\Omega$cm$^2$','Interpreter','latex')
xlabel('frequency / Hz','Interpreter','latex')
xlim([min(freqs),max(freqs)])
daspect([max(freqs),max(R),max(R)])

% === additional plots ===
% Users can create their own plots here. Any variables other than `X`, `R`,
% and `freqs` will need to be extracted from `sol`.




end