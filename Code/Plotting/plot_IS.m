function plot_IS(sol)
% A function to plot data from impedance spectroscopy simulations. The
% parameters of the Nyquist and Bode plots can easily be edited and new
% plots can be added below. `sol` is a solution structure from an IS
% simulation. This function works for full solutions or reduced solutions
% obtained using the `reduced_output` parameter. For analysing and plotting
% other sinudoidally varying outputs, see the example in FourierFit.m
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
xlim([min(freqs),max(freqs)])
ylim([min([R;0]), ceil(max(R/5))*5])
zlim([floor(min(X/5))*5 max([X;0])])
plot3(ax.XLim(2)*ones(size(R)),R,X,'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto R-X plane
plot3(freqs,ax.YLim(1)*ones(size(R)),X,'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto f-X plane
plot3(freqs,R,ax.ZLim(2)*ones(size(R)),'Color',0.6*[1,1,1],'LineWidth',1.5) % projection onto f-R plane
patch([min(freqs)*1e-3 max(freqs)*1e3 max(freqs)*1e3 min(freqs)*1e-3],min(R)*[1,1,1,1],[ax.ZLim(2),ax.ZLim(2),ax.ZLim(1),ax.ZLim(1)],...
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
daspect([max(freqs),max(R),max(R)])

% === additional plots ===
% Users can create their own plots here. Any variables other than `X`, `R`,
% and `freqs` will need to be extracted from `sol`.

% ===================== example of periodic analysis ======================
% To analyse the 'impedance' of any periodic variable in the solution
% structure, use and adapt the following code. Ensure `reduced_output` is
% set to `false` in the parameters file.
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
%     %(some quantities may require an extra phase offset of pi)
%     Vp = sol(1).impedance_protocol{5}; % get voltage amplitude
%     Z(j) = Vp/fit.Sp*exp(-i*theta_S(j)); % create 'impedance'
% end
% 
% figure()    % Nyquist plot
% plot(real(Z),-imag(Z),'-sr')
% grid on
% set(gca,'DataAspectRatio',[1 1 1]) % comment if poor view
% ylabel('-imag($V/S$)','Interpreter','latex');
% xlabel('real($V/S$)','Interpreter','latex')
% 
% figure()    % Real and imaginary frequency plot
% T = tiledlayout(2,1);
% ax1 = nexttile;
% semilogx(ax1,f,-imag(Z),'x-b')
% ax2 = nexttile;
% semilogx(ax2,f,real(Z),'x-b')
% ylabel(ax1,'-imag($V/S$)','Interpreter','latex');
% ylabel(ax2,'real($V/S$)','Interpreter','latex');
% xlabel(ax2,'frequency / Hz','Interpreter','latex')
% 
% figure()    % Amplitude and phase frequency plot
% T = tiledlayout(2,1);
% ax1 = nexttile;
% semilogx(ax1,f,Sp,'x-b')
% ax2 = nexttile;
% semilogx(ax2,f,theta_S,'x-b')
% ylabel(ax1,'amplitude $S_p$','Interpreter','latex')
% ylabel(ax2,'phase $\theta_S$','Interpreter','latex')
% xlabel(ax2,'frequency / Hz')

% =========================================================================




end