 function animate_sections(sol,sections,vidlength,filename,NameValueArgs)
% A function to animate sections of a solution from IonMonger. `sol` is a
% solution structure. `sections` is a vector of integers specifying which
% sections of the protocol should be included. `length` is the length of
% the final video in seconds. `filename` is a string that will be used to
% name the final file when saving. The optional Name-Value pair arguments
% `FrameRate` and `Size` can be used to specify the frame rate in fps and
% the size of the image in pixels. Defaults are 30fps and 1920x1080p
% resolution. If the Image Processing Toolbox is installed, the video will
% automatically play in the built-in player. Otherwise, an external player
% will be required. If the Parallel Computing Toolbox is installed, frames
% will be rendered on parallel workers.

% For example, to make a 5 second video of the reverse and forward sweeps
% after a preconditioning stage and save as 'sweep_animation' in the folder
% 'Videos', type `animate_sections(sol,[2,3],5,'Videos\sweep_animation')`
% into the command window.

% check optional Name-Value arguments
arguments
    sol struct
    sections double
    vidlength double
    filename string
    NameValueArgs.Size (1,2) {mustBeNumeric} = [1920 1080];
    NameValueArgs.FrameRate (1,1) {mustBeNumeric} = 30;
    NameValueArgs.SupressOutput = false;
end

% check sol structure
if size(sol,2)>1
    % recieved structure array from IS simulation
    error(['animate_sections was given a solution structure array from an ' ...
    'impedance spectroscopy simulation. To use animate_sections for the n-th '...
    'sample frequency solution, use `animate_sections(sol(n),...)`'])
elseif isfield(sol,'X')
    % recieved reduced solution structure from IS simulation
    error(['animate_sections was given a reduced solution structure from an ' ...
    'impedance spectroscopy simulation. To use animate_sections with an IS ' ...
    'solution, ensure reduced_output=false'])
end

tic;

if isempty(sections)
    % if sections is empty, all sections will be included
    sections = [1:(length(sol.time)-1)/100];
end
    
% create video
try
    vid = VideoWriter(filename,'MPEG-4');
catch me
    if strcmp(me.identifier,'MATLAB:audiovideo:VideoWriter:fileNotWritable')
        error(['A video file with the same name is currently open in MATLAB''s '...
            'movie player. Close this window before overwriting the video file.'])
    else
        rethrow(me)
    end
end
vid.FrameRate=NameValueArgs.FrameRate;
vid.Quality=100;

N = round(vidlength*vid.FrameRate); % number of frames

open(vid);

% indices of the start and end points
ind = (100*(sections(1)-1)+1):(100*(sections(end))+1);

% construct vector of frame times, linearly spaced
time = linspace(sol.time(ind(1)),sol.time(ind(end)),N);

% unpack necessary variables
[xE,x,xH] = struct2array(sol.vectors,{'xE','x','xH'});
[dE,dH,kE,kH,N0,G0,gamma,ni2,tor,tor3,Cn,Cp,brate,Tion,b,q,gammaE,torE,torE3,...
    gammaH,torH,torH3,brateE,brateH] = struct2array(...
    sol.params,{'dE','dH','kE','kH','N0','G0','gamma','ni2','tor','tor3',...
    'Cn','Cp','brate','Tion','b','q','gammaE','torE','torE3','gammaH', ...
    'torH','torH3','brateE','brateH'});

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

% calculate dimensional recombination rates
R = sol.params.R(n/(dE*kE),p/(dH*kH),P/N0)*G0;
Rbim = (brate*(n/(dE*kE).*p/(dH*kH)-ni2))*G0;
SRH = sol.params.SRH(n/(dE*kE),p/(dH*kH),gamma,ni2,tor,tor3)*G0;
Auger = sol.params.Auger(n/(dE*kE),p/(dH*kH),Cn,Cp,ni2)*G0;
[X,T] = meshgrid(x*1e-9/b,time/Tion); % nondimensional time and space mesh
G = sol.params.G(X,T)*G0;

% surface recombination rates are divided by perovskite width to be 
% dimensionally consistent with bulk recombination rates
Rl_SRH = sol.params.SRH(nE(:,end)/dE,p(:,1)/(dH*kH),gammaE,ni2/kE,torE,torE3)*G0;
Rr_SRH = sol.params.SRH(pH(:,1)/dH,n(:,end)/(dE*kE),gammaH,ni2/kH,torH,torH3)*G0;

Rl_bim = brateE*(nE(:,end)/dE.*p(:,1)/(dH*kH)-ni2/kE)*G0;
Rr_bim = brateH*(n(:,end)/(dE*kE).*pH(:,1)/(dH)-ni2/kH)*G0;

Rl = Rl_SRH+Rl_bim;
Rr = Rr_SRH+Rr_bim;

% create structure to contain frames
F(N) = struct('cdata',[],'colormap',[]);

% get positive elements of G and R to set correct scales for dark
% simulations or slightly negative rates
GRpos = [G,R,Rl_bim,Rl_SRH,Rr_bim,Rr_SRH]; 
GRpos(GRpos<=0) = nan;

carriers = [P,n,nE,p,pH];

% predetermine Y limits for plotting
YLims = [-2, ceil(max(J)/5)*5;
    10.^(floor(min(log10(carriers),[],'all')/2)*2),10.^(ceil(max(log10(carriers),[],'all')/2)*2);
    min([phiE,phi,phiH],[],'all')-0.1, max([phiE,phi,phiH],[],'all');
    10.^(floor(min(log10(GRpos),[],'all')/2)*2), 10.^(ceil(max(log10(GRpos),[],'all')/2)*2)];

playbackspeed = (time(end)-time(1))/(vidlength);

% package up all the data necessary to create a frame
dstrbns = struct('P',P,'phi',phi,'n',n,'p',p,'phiE',phiE,'nE',nE,'phiH',phiH,...
    'pH',pH,'R',R,'Rbim',Rbim,'SRH',SRH,'Auger',Auger,'G',G,'Rl_SRH',Rl_SRH, ...
    'Rl_bim',Rl_bim,'Rr_SRH',Rr_SRH,'Rr_bim',Rr_bim);
framedata = struct('sol',sol,'Size',NameValueArgs.Size,'dstrbns',dstrbns, ...
    'J',J,'V',V,'Rl',Rl,'Rr',Rr,'sections',sections,'time',time,'YLims', ...
    YLims,'playbackspeed',playbackspeed);

% make title frames
fig = figure('Position',[0 0 NameValueArgs.Size],'Visible','off');
ax = axes('Position',[0 0 1 1],'Color','w','Units','normalized');
imshow('Code/Plotting/IM_logo.png','Parent',ax,'Reduce',false)
frame = getframe(fig);
titleframetime = 2; % time for title frame in seconds
for i = 1:round(titleframetime*NameValueArgs.FrameRate)
    writeVideo(vid,frame)
end

frame = create_initial_frame(framedata);
for i = 1:round(titleframetime*NameValueArgs.FrameRate)
    writeVideo(vid,frame)
end

% create each frame and save to frames structure
if ~isempty(ver('parallel')) % check for parallel computing toolbox
    % parallel computing toolbox installed
    pool = gcp;
    fprintf('\nParallel computing toolbox detected \nRendering frames on %s workers \n\n', num2str(pool.NumWorkers))
    parfor i = 1:N
        F(i) = create_frame(framedata,i);
        if ~NameValueArgs.SupressOutput
            fprintf('frame %s of %s \n', num2str(i),num2str(N))
        end
    end
else
    % parallel computing toolbox not installed
    fprintf('\nParallel computing toolbox not detected \nRendering frames without parallel computing \n\n')
    for i = 1:N
        F(i) = create_frame(framedata,i);
        if ~NameValueArgs.SupressOutput
            fprintf('frame %s of %s \n', num2str(i),num2str(N))
        end
    end
end

% write each frame to video
for i = 1:N
    writeVideo(vid,F(i));
end

close(vid);   

fprintf('animation completed at %s, rendering %s frames in %ss \nvideo saved as %s.mp4 \n', ...
    datestr(now),num2str(N),num2str(toc), filename)

if ~isempty(ver('images')) % check for image processing toolbox
    str = append(filename,'.mp4');
    vidhandle = implay(str); % play video file within MATLAB
    play(vidhandle.DataSource.Controls)
    vidhandle.Parent.WindowState = 'maximized'; % play in full-screen
else
    warning('Image processing toolbox not installed. Video will need to be opened from an external player')
end

end

function frame = create_frame(framedata,i)
    % unpack data
    [sections,time,J,V,YLims,sol,Rl,Rr,Size,playbackspeed] = struct2array(framedata,{'sections','time','J','V','YLims','sol','Rl', ...
        'Rr','Size','playbackspeed'});
    [xE,x,xH] = struct2array(sol.vectors,{'xE','x','xH'});
    [P,phi,n,p,phiE,nE,phiH,pH,R,Rbim,SRH,Auger,G,Rl_SRH,Rl_bim,Rr_SRH, ...
        Rr_bim] = struct2array(framedata.dstrbns,{'P','phi','n','p','phiE', ...
        'nE','phiH','pH','R','Rbim','SRH','Auger','G','Rl_SRH','Rl_bim', ...
        'Rr_SRH','Rr_bim'});
    
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
    colororder(ax1,lines)
    start = 100*(sections(1)-1)+1;
    for j = 1:size(sections,2) % for each section
        finish = 100*sections(j)+1; % end index of section
        
        % find indices of time vector within this section
        k = find(time>=sol.time(start)& time<=sol.time(finish));
        
        % plot each current section in a different color until reaching the
        % present time
        plot(ax1,V(k(1):min([i,k(end)+1])),J(k(1):min([i,k(end)+1])),'LineWidth',1.5)
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
    plot(ax3,[xE;x;xH],[phiE(i,:),phi(i,:),phiH(i,:)],'-','LineWidth',1.5,...
        'Color',[110, 19, 128]/255)
    
    % recombination plot
    patch(ax4,x([1,end,end,1]),YLims(4,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none','HandleVisibility','off')
    if length(find([any(Rbim>0,'all') any(SRH>0,'all') any(Auger>0,'all')]))>1 % if more than one type of bulk recombination is present
        plot(ax4,x,R(i,:),'DisplayName','total bulk recombination',...
        'LineWidth',3,'Color',[12,0,120]/256)
    end
    if any(Rbim>0); plot(ax4,x,Rbim(i,:),'m','DisplayName','Bimolecular','LineWidth',2); end
    if any(SRH>0); plot(ax4,x,SRH(i,:),'-','DisplayName','SRH','LineWidth',2,'Color',[255,153,0]/256); end
    if any(Auger>0); plot(ax4,x,Auger(i,:),'-','DisplayName','Auger','LineWidth',2,'Color',[255,25,80]/256); end
    plot(ax4,x,G(i,:),'-','DisplayName','generation rate','LineWidth',2,'Color',[0,120,12]/256)
    if any(Rl_SRH>0); plot(ax4,x(1),Rl_SRH(i),'xb','DisplayName','left surface SRH','LineWidth',1.5,'MarkerSize',10); end
    if any(Rl_bim>0); plot(ax4,x(1),Rl_bim(i),'ob','DisplayName','left surface bim','LineWidth',1.5,'MarkerSize',10); end
    if any(Rr_SRH>0); plot(ax4,x(end),Rr_SRH(i),'xr','DisplayName','right surface SRH','LineWidth',1.5,'MarkerSize',10); end
    if any(Rr_bim>0); plot(ax4,x(end),Rr_bim(i),'or','DisplayName','right surface bim','LineWidth',1.5,'MarkerSize',10); end
    
    % axis labels
    xlabel(ax1,'applied voltage (V)')
    ylabel(ax1,'current density (mAcm$^{-2}$)')
    xlabel(ax2,'$x$ (nm)')
    ylabel(ax2,'number density (m$^{-3}$)')
    xlabel(ax3,'$x$ (nm)')
    ylabel(ax3,'electric potential (V)')
    xlabel(ax4,'$x$ (nm)')
    ylabel(ax4,'recombination rate (m$^{-3}$s$^{-1}$)')
    
    % add timestamp
    if log10(time(end)-time(1))<-1 % if time span is small use more precise timestamp
        prec = ceil(-log10(time(end)-time(1)))+3; % retain three significant figures
        str = ['$t =$ ' num2str(time(i), prec) 's'];
    else
        str = ['$t =$ ' datestr(seconds(time(i)), 'HH:MM:SS.FFF')];
    end
    title(T,str,'FontSize',19,'Interpreter','latex')
    txt = annotation(fig,'textbox', [0.8, 0.9, 0.1, 0.1],...
        'String', ['playing at $\times$' num2str(playbackspeed,3) ' speed'],...
        'Interpreter','latex',...
        'FontSize',18,...
        'FitBoxToText','on',...
        'EdgeColor','none',...
        'HorizontalAlignment','right');
    txt.Position(1) = 0.95 - txt.Position(3);
    
    % legends
    leg2 = legend(ax2,'NumColumns',2,'Location','south');
    leg2.BoxFace.ColorType='truecoloralpha';
    leg2.BoxFace.ColorData=uint8(255*[1 1 1 0.5]'); % set transparancy
    
    % position legend in centre of perovskite layer
    f = (-xE(1)+x(end)/2)/(-xE(1)+xH(end));
    leg2.Position(1) = ax2.Position(1) + f*ax2.Position(3)-leg2.Position(3)/2;
    
    leg4 = legend(ax4,'NumColumns',2,'Location','south');
    leg4.BoxFace.ColorType='truecoloralpha';
    leg4.BoxFace.ColorData=uint8(255*[1 1 1 0.5]'); % set transparancy
    
    % set axis properties
    set(T, 'Padding', 'compact', 'TileSpacing', 'compact')
    set(ax1,'YLim',YLims(1,:),'XLim',[floor(min(V)/0.1)*0.1,...
        ceil(max(V)/0.1)*0.1],'GridAlpha',0.4)
    set(ax2,'YScale','log','YLim',YLims(2,:),'XLim',[xE(1) xH(end)])
    set(ax3,'YLim',YLims(3,:),'XLim',[xE(1) xH(end)])
    set(ax4,'YScale','log','YLim',YLims(4,:),'XLim',[x(1) x(end)])
    
    frame = getframe(fig);
    
    close(fignum)

end

function frame = create_initial_frame(framedata)
    % unpack data
    [sections,time,J,V,YLims,sol,Rl,Rr,Size] = struct2array(framedata,{'sections','time','J','V','YLims','sol','Rl', ...
        'Rr','Size'});
    [xE,x,xH] = struct2array(sol.vectors,{'xE','x','xH'});
    [P,phi,n,p,phiE,nE,phiH,pH,R,Rbim,SRH,Auger,G,Rl_SRH,Rl_bim,Rr_SRH, ...
        Rr_bim] = struct2array(framedata.dstrbns,{'P','phi','n','p','phiE', ...
        'nE','phiH','pH','R','Rbim','SRH','Auger','G','Rl_SRH','Rl_bim', ...
        'Rr_SRH','Rr_bim'});
    
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
    
    % carrier densities plot
    patch(ax2,xE([1,end,end,1]),YLims(2,[1,1,2,2]),'b','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','ETL')
    patch(ax2,x([1,end,end,1]),YLims(2,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','Perovskite')
    patch(ax2,xH([1,end,end,1]),YLims(2,[1,1,2,2]),'r','FaceAlpha',0.15,...
        'EdgeColor','none','DisplayName','HTL')
    plot(ax2,nan,nan,'-b','DisplayName','electrons',...
        'LineWidth',1.5)
    plot(ax2,nan,nan,'DisplayName','anion vacancies','Color',[153, 125, 32]/255,...
        'LineWidth',1.5)
    plot(ax2,nan,nan,'-r','DisplayName','holes',...
        'LineWidth',1.5)
    
    % phi plot
    patch(ax3,xE([1,end,end,1]),YLims(3,[1,1,2,2]),'b','FaceAlpha',0.15,...
        'EdgeColor','none')
    patch(ax3,x([1,end,end,1]),YLims(3,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none')
    patch(ax3,xH([1,end,end,1]),YLims(3,[1,1,2,2]),'r','FaceAlpha',0.15,...
        'EdgeColor','none')
    plot(ax3,nan,nan,'-','LineWidth',1.5,...
        'Color',[110, 19, 128]/255)
    
    % recombination plot
    patch(ax4,x([1,end,end,1]),YLims(4,[1,1,2,2]),'y','FaceAlpha',0.15,...
        'EdgeColor','none','HandleVisibility','off')
    if length(find([any(Rbim>0,'all') any(SRH>0,'all') any(Auger>0,'all')]))>1 % if more than one type of bulk recombination is present
        plot(ax4,nan,nan,'DisplayName','total bulk recombination',...
        'LineWidth',3,'Color',[12,0,120]/256)
    end
    if any(Rbim>0); plot(ax4,nan,nan,'m','DisplayName','Bimolecular','LineWidth',2); end
    if any(SRH>0); plot(ax4,nan,nan,'-','DisplayName','SRH','LineWidth',2,'Color',[255,153,0]/256); end
    if any(Auger>0); plot(ax4,nan,nan,'-','DisplayName','Auger','LineWidth',2,'Color',[255,25,80]/256); end
    plot(ax4,nan,nan,'-','DisplayName','generation rate','LineWidth',2,'Color',[0,120,12]/256)
    if any(Rl_SRH>0); plot(ax4,nan,nan,'xb','DisplayName','left surface SRH','LineWidth',1.5,'MarkerSize',10); end
    if any(Rl_bim>0); plot(ax4,nan,nan,'ob','DisplayName','left surface bim','LineWidth',1.5,'MarkerSize',10); end
    if any(Rr_SRH>0); plot(ax4,nan,nan,'xr','DisplayName','right surface SRH','LineWidth',1.5,'MarkerSize',10); end
    if any(Rr_bim>0); plot(ax4,nan,nan,'or','DisplayName','right surface bim','LineWidth',1.5,'MarkerSize',10); end
    
    % axis labels
    xlabel(ax1,'applied voltage (V)')
    ylabel(ax1,'current density (mAcm$^{-2}$)')
    xlabel(ax2,'$x$ (nm)')
    ylabel(ax2,'number density (m$^{-3}$)')
    xlabel(ax3,'$x$ (nm)')
    ylabel(ax3,'electric potential (V)')
    xlabel(ax4,'$x$ (nm)')
    ylabel(ax4,'recombination rate (m$^{-3}$s$^{-1}$)')
    
    % add timestamp
    if log10(time(end)-time(1))<-1 % if time span is small use more precise timestamp
        prec = ceil(-log10(time(end)-time(1)))+3; % retain three significant figures
        str = ['$t =$ ' num2str(time(1), prec) 's'];
    else
        str = ['$t =$ ' datestr(seconds(time(1)), 'HH:MM:SS.FFF')];
    end
    title(T,str,'FontSize',19,'Interpreter','latex')
    
    % legends
    leg2 = legend(ax2,'NumColumns',2,'Location','south');
    leg2.BoxFace.ColorType='truecoloralpha';
    leg2.BoxFace.ColorData=uint8(255*[1 1 1 0.5]'); % set transparancy
    
    % position legend in centre of perovskite layer
    f = (-xE(1)+x(end)/2)/(-xE(1)+xH(end));
    leg2.Position(1) = ax2.Position(1) + f*ax2.Position(3)-leg2.Position(3)/2;
    
    leg4 = legend(ax4,'NumColumns',2,'Location','south');
    leg4.BoxFace.ColorType='truecoloralpha';
    leg4.BoxFace.ColorData=uint8(255*[1 1 1 0.5]'); % set transparancy
    
    % set axis properties
    set(T, 'Padding', 'compact', 'TileSpacing', 'compact')
    set(ax1,'YLim',YLims(1,:),'XLim',[floor(min(V)/0.1)*0.1,...
        ceil(max(V)/0.1)*0.1],'GridAlpha',0.4)
    set(ax2,'YScale','log','YLim',YLims(2,:),'XLim',[xE(1) xH(end)])
    set(ax3,'YLim',YLims(3,:),'XLim',[xE(1) xH(end)])
    set(ax4,'YScale','log','YLim',YLims(4,:),'XLim',[x(1) x(end)])
    
    text(ax1,0.5,0.5,'J-V curve','Units','normalized','FontSize',30,'HorizontalAlignment','center')
    text(ax2,0.5,0.5,'carrier densities','Units','normalized','FontSize',30,'HorizontalAlignment','center')
    text(ax3,0.5,0.5,'electric potential','Units','normalized','FontSize',30,'HorizontalAlignment','center')
    text(ax4,0.5,0.5,'recombination rates','Units','normalized','FontSize',30,'HorizontalAlignment','center')
    
    frame = getframe(fig);
    
    close(fignum)

end


