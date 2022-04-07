function plot_sections(sol,sections)
% Plots a J-V curve(s), assuming psi(t) was created using construct_potential.
% The inputs are a solution structure and an array containing a list of
% integers indicating which of the sections of the voltage protocol
% correspond to the reverse and forward scans of a J-V curve (e.g. [2,3]).
% See GUIDE.md for notes on constructing a suitable simulation protocol.
% If there is an odd number of sections, the first is assumed to correspond
% to a preconditioning step.

% Parameter input
[J, V, Jl, Jr] = struct2array(sol,{'J','V','Jl','Jr'});

% Set default figure options
set(0,'defaultAxesFontSize',18); % Make axes labels larger
set(0,'defaultTextInterpreter','latex') % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex') % For latex tick labels
set(0,'defaultLegendInterpreter','latex') % For latex legends

% Figure setup
figure;
hold on;
ylim([-15, 30]);
ylabel('Current density (mA$\cdot$cm$^{-2}$)');
xlim([0, 1.2]);
xlabel('Voltage (V)');
line([0,1.2],[0,0],'Color','k','HandleVisibility','off'); % x-axis, i.e. J=0
set(gca,'Layer','top');
grid on;
box on;

% Clever colour picking
colour{1} = [0.6, 0, 0.9]; % purple
colour{2} = [0, 0.9, 0.6]; % green
colour{3} = [0.9, 0.6, 0]; % orange

% Plot each dataset
odd = mod(length(sections),2);
if odd % odd number of sections
    % Plot the first section
    start = (sections(1)-1)*100+1;
    plot(V(start:start+100),J(start:start+100),'-','Color',[0, 0, 0]);
    plot(V(start:start+100),Jl(start:start+100),':b','HandleVisibility','off');
    plot(V(start:start+100),Jr(start:start+100),':r','HandleVisibility','off');
end
for i = 1:floor(length(sections)/2)
    % Plot first scan (usually reverse)
    start = (sections(odd+2*i-1)-1)*100+1;
    plot(V(start:start+100),J(start:start+100),'-','Color',colour{i});
    plot(V(start:start+100),Jl(start:start+100),'-.b','HandleVisibility','off');
    plot(V(start:start+100),Jr(start:start+100),'-.r','HandleVisibility','off');
    % Plot second scan (usually forward)
    start = (sections(odd+2*i)-1)*100+1;
    plot(V(start:start+100),J(start:start+100),'--','Color',colour{i}); %,'HandleVisibility','off');
    plot(V(start:start+100),Jl(start:start+100),':b','HandleVisibility','off');
    plot(V(start:start+100),Jr(start:start+100),':r','HandleVisibility','off');
end

% Add a legend
legnames = {sprintf('Section %g',sections(1))};
for i = 2:length(sections)
    legnames = {legnames{:}, sprintf('Section %g',sections(i))};
end
legend(legnames,'Location','Best','FontSize',12);

end
