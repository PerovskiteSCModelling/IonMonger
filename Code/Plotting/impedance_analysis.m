function [X,R] = impedance_analysis(sol)

nf = sol(1).impedance_protocol{7}; % number of frequencies to be sampled
min_f = sol(1).impedance_protocol{2};
V0 = sol(1).impedance_protocol{4};
max_f = sol(1).impedance_protocol{3};
n_wave = sol(1).impedance_protocol{8};

freqs = logspace(log10(min_f),log10(max_f),nf);

Z = nan(nf,1)+i*nan(nf,1);
for j = 1:nf
    try
        Z(j) = Z_from_J(sol(j));
    end
end

R = real(Z);
X = imag(Z);

% make plots

figure(1)
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
title(ax1,['$V_{DC} =$ ' num2str(V0) 'V'],'Interpreter','latex')

set(gcf,'Units','pixels','Position',[200,200,1300,600])

T.TileSpacing = 'compact';
drawnow;

% Nyquist plot

figure(2)
plot(R,-X,'-sr')
grid on
set(gca,'TickLabelInterpreter','latex','FontSize',18)
ylabel('-X / $\Omega$cm$^2$','Interpreter','latex')
xlabel('R / $\Omega$cm$^2$','Interpreter','latex')
title(['$V_{DC} =$ ' num2str(V0) 'V'],'Interpreter','latex')
drawnow;

end
