function animate_sections(sol,sections,length,filename,NameValueArgs)
% A function to animate sections of a solution from IonMonger. `sol` is a
% solution structure. `sections` is a vector of integers specifying which
% sections of the protocol should be included. `length` is the length of
% the final video in seconds. `filename` is a string that will be used to
% name the final file when saving. The optional Name-Value pair arguments
% `FrameRate` and `Size` can be used to specify the frame rate in fps and
% the size of the image in pixels. Defaults are 30fps and 1920x1080p
% resolution.

% For example, to make a 5 second video of the reverse and forward sweeps
% after a preconditioning stage and save as 'sweep_animation' in the folder
% 'Videos', type `animate_sections(sol,[2:3],5,'Videos\sweep_animation')`
% into the command window.

% check optional Name-Value arguments
arguments
    sol struct
    sections double
    length double
    filename string
    NameValueArgs.Size (1,2) {mustBeNumeric} = [1920 1080];
    NameValueArgs.FrameRate (1,1) {mustBeNumeric} = 30;
end

tic;

% create video
vid = VideoWriter(filename,'MPEG-4');
vid.FrameRate=NameValueArgs.FrameRate;
vid.Quality=100;

N = round(length*vid.FrameRate); % number of frames

open(vid);

% indices of the start and end points
ind = (100*(sections(1)-1)+1):(100*(sections(end))+1);

% construct vector of frame times
time = linspace(sol.time(ind(1)),sol.time(ind(end)),N);

% unpack necessary variables
[xE,x,xH] = struct2array(sol.vectors,{'xE','x','xH'});
[dE,dH,kE,kH,N0,G0,gamma,ni2,tor,tor3,Cn,Cp,brate,Tion,b,q] = struct2array(...
    sol.params,{'dE','dH','kE','kH','N0','G0','gamma','ni2','tor','tor3',...
    'Cn','Cp','brate','Tion','b','q'});

% interpolate solutions onto time grid
P = interp1(sol.time,sol.dstrbns.P,time);
phi = interp1(sol.time,sol.dstrbns.phi,time);
n = interp1(sol.time,sol.dstrbns.n,time);
p = interp1(sol.time,sol.dstrbns.p,time);
nE = interp1(sol.time,sol.dstrbns.nE,time);
phiE = interp1(sol.time,sol.dstrbns.phiE,time);
pH = interp1(sol.time,sol.dstrbns.pH,time);
phiH = interp1(sol.time,sol.dstrbns.phiH,time);

J = interp1(sol.time,sol.J,time);
V = interp1(sol.time,sol.V,time);
Rl = interp1(sol.time,-sol.Jl/(b*q),time);
Rr = interp1(sol.time,-sol.Jr/(b*q),time);

R = sol.params.R(n/(dE*kE),p/(dH*kH),P/N0)*G0;
Rbim = (brate*(n/(dE*kE).*p/(dH*kH)-ni2))*G0;
SRH = sol.params.SRH(n/(dE*kE),p/(dH*kH),gamma,ni2,tor,tor3)*G0;
Auger = sol.params.Auger(n/(dE*kE),p/(dH*kH),Cn,Cp,ni2)*G0;
[X,T] = meshgrid(x*1e-9/b,time/Tion); % nondimensional time and space mesh
G = sol.params.G(X,T)*G0;

% create structure to contain frames
F(N) = struct('cdata',[],'colormap',[]);

% predetermine Y limits for plotting
YLims = [-2, ceil(max(J)/5)*5;
            0.1*min([P,n,nE,p,pH],[],'all'),10*max([P,n,nE,p,pH],[],'all');
            min([phiE,phi,phiH],[],'all')-0.1, max([phiE,phi,phiH],[],'all');
            min([R,G],[],'all'), 10.^(ceil(max(log10([R,G]),[],'all')))];
        
% package up all the data necessary to create a frame
dstrbns = struct('P',P,'phi',phi,'n',n,'p',p,'phiE',phiE,'nE',nE,'phiH',phiH,...
    'pH',pH,'R',R,'Rbim',Rbim,'SRH',SRH,'Auger',Auger,'G',G);
framedata = struct('sol',sol,'Size',NameValueArgs.Size,'dstrbns',dstrbns,'J',J,'V',V,'Rl',Rl,'Rr',Rr,'sections',...
    sections,'time',time,'YLims',YLims);

% create each frame and save to frames structure
if ~isempty(ver('parallel')) % check for parallel computing toolbox
    % parallel computing toolbox installed
    pool = gcp;
    fprintf('\nParallel computing toolbox detected \nRendering frames on %s workers \n\n', num2str(pool.NumWorkers))
    parfor i = 1:N
        F(i) = create_frame(framedata,i);
        fprintf('frame %s of %s \n', num2str(i),num2str(N))
    end
else
    % parallel computing toolbox not installed
    fprintf('\nParallel computing toolbox not detected \nRendering frames without parallel computing \n\n')
    for i = 1:N
        F(i) = create_frame(framedata,i);
        fprintf('frame %s of %s \n', num2str(i),num2str(N))
    end
end

% write each frame to video
for i = 1:N
    writeVideo(vid,F(i));
end

close(vid);   

fprintf('animation completed at %s, rendering %s frames in %ss \nvideo saved as %s \n', ...
    datestr(now),num2str(N),num2str(toc), filename)

end

function frame = create_frame(framedata,i)

    % unpack data
    [sections,time,J,V,YLims,sol,Rl,Rr,Size] = struct2array(framedata,{'sections','time','J','V','YLims','sol','Rl','Rr','Size'});
    [xE,x,xH] = struct2array(sol.vectors,{'xE','x','xH'});
    [P,phi,n,p,phiE,nE,phiH,pH,R,Rbim,SRH,Auger,G] = struct2array(framedata.dstrbns,...
        {'P','phi','n','p','phiE','nE','phiH','pH','R','Rbim','SRH','Auger','G'});
    
    % Set default figure options
    set(0,'defaultAxesFontSize',18); % Make axes labels larger
    set(0,'defaultTextInterpreter','latex') % For latex axis labels
    set(0,'defaultAxesTickLabelInterpreter','latex') % For latex tick labels
    set(0,'defaultLegendInterpreter','latex') % For latex legends

    fig = figure('Position',[0 0 Size(1) Size(2)],'Visible','off');
    fignum = fig.Number;
    clf(fignum)
    T=tiledlayout(2,2);
    ax1=nexttile(1); % jv
    ax2=nexttile(2); % densities
    ax3=nexttile(3); % potential
    ax4=nexttile(4); % recombination
    
    hold(ax1,'on')
    hold(ax2,'on')
    hold(ax3,'on')
    hold(ax4,'on')
    box(ax1,'on')
    box(ax2,'on')
    box(ax3,'on')
    box(ax4,'on')
    grid(ax1,'on')
    grid(ax3,'on')
    
    % current plot
    start = 100*(sections(1)-1)+1;
    for j = 1:size(sections,2) % for each section
        finish = 100*sections(j)+1; % end index of section
    
        % find indices of time vector within this section
        k = find(time>=sol.time(start)& time<=sol.time(finish));
        
        % plot each current section in a different color until reaching the
        % present time
        plot(ax1,V(k(1):min([i,k(end)])),J(k(1):min([i,k(end)])),'LineWidth',1.5)
        plot(ax1,V(i)*ones(1,2),YLims(1,:),'-k')
        
        start = finish; % last index becomes the first index of next line
    end
    
    % carrier densities plot
    patch(ax2,xE([1,end,end,1]),YLims(2,[1,1,2,2]),'b','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','ETL')
    patch(ax2,x([1,end,end,1]),YLims(2,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','Perovskite')
    patch(ax2,xH([1,end,end,1]),YLims(2,[1,1,2,2]),'r','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','HTL')
    plot(ax2,[xE;,x],[nE(i,:),n(i,:)],'-b','DisplayName','electrons',...
        'LineWidth',1.5)
    plot(ax2,x,P(i,:),'DisplayName','anion vacancies','Color',[153, 125, 32]/255,...
        'LineWidth',1.5)
    plot(ax2,[x; xH],[p(i,:),pH(i,:)],'-r','DisplayName','holes',...
        'LineWidth',1.5)
    
    % phi plot
    patch(ax3,xE([1,end,end,1]),YLims(3,[1,1,2,2]),'b','FaceAlpha',0.15,...
        'EdgeColor','none')
    patch(ax3,x([1,end,end,1]),YLims(3,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none')
    patch(ax3,xH([1,end,end,1]),YLims(3,[1,1,2,2]),'r','FaceAlpha',0.15,...
        'EdgeColor','none')
    plot(ax3,[xE;x;xH],[phiE(i,:),phi(i,:),phiH(i,:)],'m','LineWidth',1.5)
    
    % recombination plot
    patch(ax4,x([1,end,end,1]),YLims(4,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none','HandleVisibility','off')
    plot(ax4,x,R(i,:),'g','DisplayName','total recombination rate',...
        'LineWidth',3)
    plot(ax4,x,Rbim(i,:),'m','DisplayName','Bimolecular','LineWidth',1)
    plot(ax4,x,SRH(i,:),'r','DisplayName','SRH','LineWidth',1)
    plot(ax4,x,Auger(i,:),'b','DisplayName','Auger','LineWidth',1)
    plot(ax4,x,G(i,:),'--k','DisplayName','generation rate','LineWidth',1.5)
    plot(ax4,x(1),Rl(i),'ob','DisplayName','left surface $\frac{R_l}{b}$','LineWidth',1.5,'MarkerSize',10)
    plot(ax4,x(end),Rr(i),'or','DisplayName','right surface $\frac{R_r}{b}$','LineWidth',1.5,'MarkerSize',10)
    
    % axis labels
    xlabel(ax1,'applied voltage (V)')
    ylabel(ax1,'current density (mAcm$^{-2}$)')
    xlabel(ax2,'$x$ (nm)')
    ylabel(ax2,'number density (m$^{-3}$)')
    xlabel(ax3,'$x$ (nm)')
    ylabel(ax3,'electric potential (V)')
    xlabel(ax4,'$x$ (nm)')
    ylabel(ax4,'recombination rate (m$^{-3}$s$^{-1}$)')
    
    title(T,['$t =$ ' datestr(seconds(time(i)), 'HH:MM:SS.FFF')],...
        'FontSize',19,'Interpreter','latex')
    
    % legends
    leg2 = legend(ax2,'NumColumns',2,'Location','south');
    % position legend in centre of perovskite layer
    f = (-xE(1)+x(end)/2)/(-xE(1)+xH(end));
    leg2.Position(1) = ax2.Position(1) + f*ax2.Position(3)-leg2.Position(3)/2;
    leg4 = legend(ax4,'NumColumns',2,'Location','south');
    
    set(T, 'Padding', 'compact', 'TileSpacing', 'compact')
    set(ax1,'YLim',YLims(1,:),'XLim',[floor(min(V)/0.1)*0.1,...
        ceil(max(V)/0.1)*0.1],'GridAlpha',0.4)
    set(ax2,'YScale','log','YLim',YLims(2,:),'XLim',[xE(1) xH(end)])
    set(ax3,'YLim',YLims(3,:),'XLim',[xE(1) xH(end)])
    set(ax4,'YScale','log','YLim',YLims(4,:),'XLim',[x(1) x(end)])
    
    frame = getframe(fig);
    
    close(fignum)

end
