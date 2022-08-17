function animate_sections(sol,sections,vidlength,filename,NameValueArgs)
% A function to animate sections of a solution from IonMonger. `sol` is a
% solution structure. `sections` is a vector of integers specifying which
% sections of the protocol should be included. `vidlength` is the length of
% the final video in seconds. `filename` is a string that will be used to
% name the final file when saving. The optional Name-Value pair arguments
% `FrameRate` and `Size` can be used to specify the frame rate in fps and
% the size of the image in pixels. Defaults are 30fps and 1280x720p
% resolution (choose 1920x1080p only if your screen has a strictly higher
% resolution than this). If the Image Processing Toolbox is installed, the
% video will automatically play in the built-in player. Otherwise, an
% external player will be required. If the Parallel Computing Toolbox is
% installed, frames will be rendered on parallel workers.

% For example, to make a 5 second video of the reverse and forward sweeps
% after a preconditioning stage and save as 'sweep_animation' in the folder
% 'Videos', type `animate_sections(sol,[2,3],5,'Videos\sweep_animation')`
% into the command window.

%% Extract data and setup

% Check optional Name-Value arguments
arguments
    sol struct
    sections double
    vidlength double
    filename string
    NameValueArgs.Size (1,2) {mustBeNumeric} = [1280 720]; % [1920 1080];
    NameValueArgs.FrameRate (1,1) {mustBeNumeric} = 30;
    NameValueArgs.Verbose = true;
end

% Check sol structure
if size(sol,2)>1 % received structure array from IS simulation
    error(['animate_sections was given a solution structure array from an ' ...
        'impedance spectroscopy simulation. To use animate_sections for the '...
        'n-th sample frequency solution, use `animate_sections(sol(n),...)`']);
elseif isfield(sol,'X') % received reduced solution structure from IS simulation
    error(['animate_sections was given a reduced solution structure from an ' ...
        'impedance spectroscopy simulation. To use animate_sections with an ' ...
        'IS solution, ensure reduced_output=false']);
end

tic;

% Unpack parameters, vectors and dimensional solution variables
[Size, FrameRate, Verbose] = ...
    struct2array(NameValueArgs, {'Size','FrameRate','Verbose'});
[n0, p0, N0, G0, gamma, ni2, tor, tor3, Cn, Cp, brate, Tion, b, gammaE, ...
    torE, torE3, gammaH, torH, torH3, brateE, brateH, G, R, SRH, Auger, ...
    light] = ...
    struct2array(sol.params, {'n0','p0','N0','G0','gamma','ni2','tor', ...
    'tor3','Cn','Cp','brate','Tion','b','gammaE','torE','torE3','gammaH', ...
    'torH','torH3','brateE','brateH','G','R','SRH','Auger','light'});
[time, J, V, vectors] = struct2array(sol, {'time','J','V','vectors'});
[x] = struct2array(vectors,{'x'});
[P, phi, n, p, phiE, nE, phiH, pH] = ...
    struct2array(sol.dstrbns,{'P','phi','n','p','phiE','nE','phiH','pH'});

% If sections is empty, include all sections
if isempty(sections)
    sections = [1:(length(time)-1)/100];
end

% Create video object
try
    vid = VideoWriter(filename,'MPEG-4');
catch me
    if strcmp(me.identifier,'MATLAB:audiovideo:VideoWriter:fileNotWritable')
        error(['A video file with the same name is currently open in MATLAB''s '...
            'movie player. Close this window before overwriting the video file.']);
    else
        rethrow(me);
    end
end
vid.FrameRate = FrameRate;
vid.Quality = 100;
Num = round(vidlength*vid.FrameRate); % number of frames

% Open video object
open(vid);

% Set indices of the start and end points
ind = (100*(sections(1)-1)+1):(100*(sections(end))+1);

% Construct vector of frame times, linearly spaced
frame_times = linspace(time(ind(1)),time(ind(end)),Num);

% Interpolate solutions onto time grid
P    = interp1(time,P,   frame_times);
phi  = interp1(time,phi, frame_times);
n    = interp1(time,n,   frame_times);
p    = interp1(time,p,   frame_times);
nE   = interp1(time,nE,  frame_times);
phiE = interp1(time,phiE,frame_times);
pH   = interp1(time,pH,  frame_times);
phiH = interp1(time,phiH,frame_times);
J    = interp1(time,J,   frame_times);
V    = interp1(time,V,   frame_times);


%% Compute recombination and generation rates

% Convert to dimensionless variables (see the scaling in numericalsolver.m)
n = n/n0;
p = p/p0;
P = P/N0;

% Calculate dimensional recombination rates (see nondimensionalise.m)
R_tot = G0*R(n,p,P);
R_bim = G0*brate*(n.*p-ni2);
R_SRH = G0*SRH(n,p,gamma,ni2,tor,tor3);
R_Aug = G0*Auger(n,p,Cn,Cp,ni2);

% Compute dimensional generation rate on a nondimensional time and space
% mesh (see the scaling in numericalsolver.m)
[X,T] = meshgrid(x*1e-9/b,frame_times/Tion);
G_tot = G(X,T)*G0;

% Divide surface recombination rates by perovskite width to be 
% dimensionally consistent with bulk recombination rates
Rl_SRH = G0*SRH(n(:,1),p(:,1),gammaE,ni2,torE,torE3);
Rr_SRH = G0*SRH(p(:,end),n(:,end),gammaH,ni2,torH,torH3);
Rl_bim = G0*brateE*(n(:,1).*p(:,1)-ni2);
Rr_bim = G0*brateH*(n(:,end).*p(:,end)-ni2);
Rl_tot = Rl_SRH+Rl_bim;
Rr_tot = Rr_SRH+Rr_bim;

% Rescale back to dimensional variables
n = n0*n;
p = p0*p;
P = N0*P;


%% Package up data and generate frames

% Create structure to contain frames
F(Num) = struct('cdata',[],'colormap',[]);

% Get positive elements of G and R to set correct scales for dark
% simulations or slightly negative rates
GRpos = [G_tot,R_tot,Rl_bim,Rl_SRH,Rr_bim,Rr_SRH]; 
GRpos(GRpos<=0) = nan;

% Compile densities into a vector
carriers = [P,n,nE,p,pH];

% Predetermine Y limits for plotting
YLims = [-2, ceil(max(J)/5)*5;
    10.^(floor(min(log10(carriers),[],'all')/2)*2),10.^(ceil(max(log10(carriers),[],'all')/2)*2);
    min([phiE,phi,phiH],[],'all')-0.1, max([phiE,phi,phiH],[],'all');
    10.^(floor(min(log10(GRpos),[],'all')/2)*2), 10.^(ceil(max(log10(GRpos),[],'all')/2)*2)];

% Set playback speed
playbackspeed = (frame_times(end)-frame_times(1))/(vidlength);

% Package up all the data necessary to create a frame
dstrbns = struct('P',P, 'phi',phi, 'n',n, 'p',p, 'phiE',phiE, 'nE',nE, ...
    'phiH',phiH, 'pH',pH, 'R_tot',R_tot, 'R_bim',R_bim, 'R_SRH',R_SRH, ...
    'R_Aug',R_Aug, 'G_tot',G_tot, 'Rl_SRH',Rl_SRH, 'Rl_bim',Rl_bim, ...
    'Rr_SRH',Rr_SRH, 'Rr_bim',Rr_bim);
framedata = struct('vectors',vectors, 'Size',Size, 'dstrbns',dstrbns, ...
    'J',J, 'V',V, 'Rl_tot',Rl_tot, 'Rr_tot',Rr_tot, 'sections',sections, ...
    'time',time, 'frame_times',frame_times, 'YLims',YLims, ...
    'playbackspeed',playbackspeed);

% Make title frames
fig = figure('Position',[0 0 Size],'Visible','off');
frame = getframe(fig);
ax = axes('Position',[0 0 1 1],'Color','w','Units','normalized');
imshow('Code/Plotting/IM_logo.png','Parent',ax,'Reduce',false);
frame = getframe(fig);
close(fig.Number);
titleframetime = 2; % time for title frame in seconds
for i = 1:round(titleframetime*FrameRate)
    writeVideo(vid,frame);
end

% Plot experimental protocol
fig = figure('Position',[0 0 Size],'Visible','off');
ax = axes('Position',[0.2 0.2 0.6 0.6],'Color','w','Units','normalized');
yyaxis left;
plot(time,sol.V,'-b','HandleVisibility','off','LineWidth',1.5);
Ys = [floor(min([sol.V; light(time')])*10)/10, ...
      ceil(max([sol.V; light(time')])*10)/10];
ylim(Ys);
set(gca,'FontSize',20,'YColor','b');
ylabel('applied voltage');
yyaxis right;
plot(time,light(time),'-r','HandleVisibility','off','LineWidth',1.5);
set(gca,'FontSize',20,'YColor','r');
ylim(Ys);
ylabel({'light intensity (Sun equiv.)'});
xlabel('time (s)');
title('Experimental protocol','Interpreter','latex');

% Create frame showing experimental protocol
patch([frame_times(1) frame_times(end) frame_times(end) frame_times(1)], ...
        [Ys([1,1,2,2])],'g','FaceAlpha',0.2,'EdgeColor','none', ...
        'DisplayName','render region');
legend('Location','best');
frame = getframe(fig);
close(fig.Number);
for i = 1:round(titleframetime*FrameRate)
    writeVideo(vid,frame);
end

% Create initial frame
frame = create_initial_frame(framedata);
for i = 1:round(titleframetime*FrameRate)
    writeVideo(vid,frame);
end

% Create each frame and save to frames structure
if ~isempty(ver('parallel')) % check for parallel computing toolbox
    % parallel computing toolbox installed
    pool = gcp;
    fprintf('\nParallel computing toolbox detected \nRendering frames on %s workers \n\n', num2str(pool.NumWorkers))
    parfor i = 1:Num
        F(i) = create_frame(framedata,i);
        if Verbose
            fprintf('frame %s of %s \n', num2str(i),num2str(Num));
        end
    end
else
    % parallel computing toolbox not installed
    fprintf('\nParallel computing toolbox not detected \nRendering frames without parallel computing \n\n')
    for i = 1:Num
        F(i) = create_frame(framedata,i);
        if Verbose
            fprintf('frame %s of %s \n', num2str(i),num2str(Num));
        end
    end
end

% Write each frame to video
for i = 1:Num
    writeVideo(vid,F(i));
end

% Complete
close(vid);
fprintf('Animation completed at %s, rendering %s frames in %ss \nvideo saved as %s.mp4 \n', ...
    datestr(now),num2str(Num),num2str(toc), filename);

% Play video file if possible
if ~isempty(ver('images')) % check for image processing toolbox
    str = append(filename,'.mp4');
    vidhandle = implay(str); % play video file within MATLAB
    play(vidhandle.DataSource.Controls);
    vidhandle.Parent.WindowState = 'maximized'; % play in full-screen
else
    warning(['Image processing toolbox not installed. Video will need ' ...
        'to be opened from an external player']);
end

end


%% Functions used by the code above:

function frame = create_frame(framedata,i)

% Unpack data
[sections, time, frame_times, J, V, YLims, vectors, Size, playbackspeed] = ...
    struct2array(framedata,{'sections','time','frame_times','J','V', ...
                            'YLims','vectors','Size','playbackspeed'});
[xE, x, xH] = struct2array(vectors,{'xE','x','xH'});
[P, phi, n, p, phiE, nE, phiH, pH, R_tot, R_bim, R_SRH, R_Aug, G_tot, ...
    Rl_SRH, Rl_bim, Rr_SRH, Rr_bim] = ...
    struct2array(framedata.dstrbns,{'P','phi','n','p','phiE','nE','phiH', ...
    'pH','R_tot','R_bim','R_SRH','R_Aug','G_tot','Rl_SRH','Rl_bim', ...
    'Rr_SRH','Rr_bim'});

% Set default figure options
set(0,'defaultAxesFontSize',18); % Make axes labels larger
set(0,'defaultTextInterpreter','latex'); % For latex axis labels
set(0,'defaultAxesTickLabelInterpreter','latex'); % For latex tick labels
set(0,'defaultLegendInterpreter','latex'); % For latex legends

% Create figure
fig = figure('Position',[0 0 Size(1) Size(2)],'Visible','off');
fignum = fig.Number;
clf(fignum);
T = tiledlayout(2,2);
ax1 = nexttile(1); % jv
ax2 = nexttile(2); % densities
ax3 = nexttile(3); % potential
ax4 = nexttile(4); % recombination
hold(ax1,'on');
hold(ax2,'on');
hold(ax3,'on');
hold(ax4,'on');
box(ax1,'on');
box(ax2,'on');
box(ax3,'on');
box(ax4,'on');
grid(ax1,'on');
grid(ax3,'on');

% Plot the current
colororder(ax1,lines);
start = 100*(sections(1)-1)+1;
for j = 1:size(sections,2) % for each section
    finish = 100*sections(j)+1; % end index of section

    % Find indices of time vector within this section
    k = find(frame_times>=time(start)& frame_times<=time(finish));

    % Plot each current section in a different color until reaching the
    % present time
    plot(ax1,V(k(1):min([i,k(end)+1])),J(k(1):min([i,k(end)+1])),'LineWidth',1.5);
    plot(ax1,V(i)*ones(1,2),YLims(1,:),'-k');

    start = finish; % last index becomes the first index of next line
end
plot(ax1,V(min([i,k(end)+1])),J(min([i,k(end)+1])),'or','LineWidth',1.5);

% Plot the carrier densities
patch(ax2,xE([1,end,end,1]),YLims(2,[1,1,2,2]),'b','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','ETL');
patch(ax2,x([1,end,end,1]),YLims(2,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','Perovskite');
patch(ax2,xH([1,end,end,1]),YLims(2,[1,1,2,2]),'r','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','HTL');
plot(ax2,[xE;x],[nE(i,:),n(i,:)],'-b','DisplayName','electrons',...
    'LineWidth',1.5);
plot(ax2,x,P(i,:),'DisplayName','anion vacancies','Color',[153, 125, 32]/255,...
    'LineWidth',1.5);
plot(ax2,[x;xH],[p(i,:),pH(i,:)],'-r','DisplayName','holes',...
    'LineWidth',1.5);

% Plot the electric potential
patch(ax3,xE([1,end,end,1]),YLims(3,[1,1,2,2]),'b','FaceAlpha',0.15,...
    'EdgeColor','none')
patch(ax3,x([1,end,end,1]),YLims(3,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none')
patch(ax3,xH([1,end,end,1]),YLims(3,[1,1,2,2]),'r','FaceAlpha',0.15,...
    'EdgeColor','none')
plot(ax3,[xE;x;xH],[phiE(i,:),phi(i,:),phiH(i,:)],'-','LineWidth',1.5,...
    'Color',[110, 19, 128]/255)

% Plot the recombination
patch(ax4,x([1,end,end,1]),YLims(4,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none','HandleVisibility','off')
if length(find([any(R_bim>0,'all') any(R_SRH>0,'all') any(R_Aug>0,'all')]))>1 % if more than one type of bulk recombination is present
    plot(ax4,x,R_tot(i,:),'DisplayName','total bulk recombination',...
    'LineWidth',3,'Color',[12,0,120]/256)
end
if any(R_bim>0); plot(ax4,x,R_bim(i,:),'m','DisplayName','Bimolecular','LineWidth',2); end
if any(R_SRH>0); plot(ax4,x,R_SRH(i,:),'-','DisplayName','SRH','LineWidth',2,'Color',[255,153,0]/256); end
if any(R_Aug>0); plot(ax4,x,R_Aug(i,:),'-','DisplayName','Auger','LineWidth',2,'Color',[255,25,80]/256); end
plot(ax4,x,G_tot(i,:),'-','DisplayName','generation rate','LineWidth',2,'Color',[0,120,12]/256)
if any(Rl_SRH>0); plot(ax4,x(1),Rl_SRH(i),'xb','DisplayName','left surface SRH','LineWidth',1.5,'MarkerSize',10); end
if any(Rl_bim>0); plot(ax4,x(1),Rl_bim(i),'ob','DisplayName','left surface bim','LineWidth',1.5,'MarkerSize',10); end
if any(Rr_SRH>0); plot(ax4,x(end),Rr_SRH(i),'xr','DisplayName','right surface SRH','LineWidth',1.5,'MarkerSize',10); end
if any(Rr_bim>0); plot(ax4,x(end),Rr_bim(i),'or','DisplayName','right surface bim','LineWidth',1.5,'MarkerSize',10); end

% Add axis labels
xlabel(ax1,'applied voltage (V)');
ylabel(ax1,'current density (mAcm$^{-2}$)');
xlabel(ax2,'$x$ (nm)');
ylabel(ax2,'number density (m$^{-3}$)');
xlabel(ax3,'$x$ (nm)');
ylabel(ax3,'electric potential (V)');
xlabel(ax4,'$x$ (nm)');
ylabel(ax4,'recombination rate (m$^{-3}$s$^{-1}$)');

% Add timestamp
if log10(frame_times(end)-frame_times(1))<-1 % if time span is small use more precise timestamp
    prec = ceil(-log10(frame_times(end)-frame_times(1)))+3; % retain three significant figures
    str = ['$t =$ ' num2str(frame_times(i), prec) 's'];
else
    str = ['$t =$ ' datestr(seconds(frame_times(i)), 'HH:MM:SS.FFF')];
end
title(T,str,'FontSize',19,'Interpreter','latex');
txt = annotation(fig,'textbox', [0.8, 0.9, 0.1, 0.1],...
    'String', ['playing at $\times$' num2str(playbackspeed,3) ' speed'],...
    'Interpreter','latex',...
    'FontSize',18,...
    'FitBoxToText','on',...
    'EdgeColor','none',...
    'HorizontalAlignment','right');
txt.Position(1) = 0.95 - txt.Position(3);

% Add legends, position legend in centre of perovskite layer
leg2 = legend(ax2,'NumColumns',2,'Location','south');
leg2.BoxFace.ColorType = 'truecoloralpha';
leg2.BoxFace.ColorData = uint8(255*[1 1 1 0.5]'); % set transparancy
f = (-xE(1)+x(end)/2)/(-xE(1)+xH(end));
leg2.Position(1) = ax2.Position(1) + f*ax2.Position(3)-leg2.Position(3)/2;
leg4 = legend(ax4,'NumColumns',2,'Location','south');
leg4.BoxFace.ColorType = 'truecoloralpha';
leg4.BoxFace.ColorData = uint8(255*[1 1 1 0.5]'); % set transparancy

% Set axis properties
set(T, 'Padding', 'compact', 'TileSpacing', 'compact');
set(ax1,'YLim',YLims(1,:),'XLim',[floor(min([0,V])/0.1)*0.1,...
    ceil(max(V)/0.1)*0.1],'GridAlpha',0.4);
set(ax2,'YScale','log','YLim',YLims(2,:),'XLim',[xE(1) xH(end)]);
set(ax3,'YLim',YLims(3,:),'XLim',[xE(1) xH(end)]);
set(ax4,'YScale','log','YLim',YLims(4,:),'XLim',[x(1) x(end)]);

% Output frame and close
frame = getframe(fig);
close(fignum);

end

function frame = create_initial_frame(framedata)

% Unpack data
[frame_times, V, YLims, vectors, Size] = ...
    struct2array(framedata,{'frame_times','V','YLims','vectors','Size'});
[xE, x, xH] = struct2array(vectors,{'xE','x','xH'});
[R_bim, R_SRH, R_Aug, Rl_SRH, Rl_bim, Rr_SRH, Rr_bim] = ...
    struct2array(framedata.dstrbns,{'R_bim','R_SRH','R_Aug','Rl_SRH', ...
                                    'Rl_bim','Rr_SRH','Rr_bim'});

% Create figure
fig = figure('Position',[0 0 Size(1) Size(2)],'Visible','off');
fignum = fig.Number;
clf(fignum);
T = tiledlayout(2,2);
ax1 = nexttile(1); % jv
ax2 = nexttile(2); % densities
ax3 = nexttile(3); % potential
ax4 = nexttile(4); % recombination
hold(ax1,'on');
hold(ax2,'on');
hold(ax3,'on');
hold(ax4,'on');
box(ax1,'on');
box(ax2,'on');
box(ax3,'on');
box(ax4,'on');

% Plot carrier densities
patch(ax2,xE([1,end,end,1]),YLims(2,[1,1,2,2]),'b','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','ETL');
patch(ax2,x([1,end,end,1]),YLims(2,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','Perovskite');
patch(ax2,xH([1,end,end,1]),YLims(2,[1,1,2,2]),'r','FaceAlpha',0.15,...
    'EdgeColor','none','DisplayName','HTL');
plot(ax2,nan,nan,'-b','DisplayName','electrons','LineWidth',1.5);
plot(ax2,nan,nan,'DisplayName','anion vacancies','LineWidth',1.5,...
    'Color',[153, 125, 32]/255);
plot(ax2,nan,nan,'-r','DisplayName','holes','LineWidth',1.5);

% Plot electric potential
patch(ax3,xE([1,end,end,1]),YLims(3,[1,1,2,2]),'b','FaceAlpha',0.15,...
    'EdgeColor','none');
patch(ax3,x([1,end,end,1]),YLims(3,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none');
patch(ax3,xH([1,end,end,1]),YLims(3,[1,1,2,2]),'r','FaceAlpha',0.15,...
    'EdgeColor','none');
plot(ax3,nan,nan,'-','LineWidth',1.5,'Color',[110, 19, 128]/255);

% Plot recombination
patch(ax4,x([1,end,end,1]),YLims(4,[1,1,2,2]),'y','FaceAlpha',0.15,...
    'EdgeColor','none','HandleVisibility','off');
if length(find([any(R_bim>0,'all') any(R_SRH>0,'all') any(R_Aug>0,'all')]))>1 % if more than one type of bulk recombination is present
    plot(ax4,nan,nan,'DisplayName','total bulk recombination',...
    'LineWidth',3,'Color',[12,0,120]/256);
end
if any(R_bim>0); plot(ax4,nan,nan,'m','DisplayName','Bimolecular','LineWidth',2); end
if any(R_SRH>0); plot(ax4,nan,nan,'-','DisplayName','SRH','LineWidth',2,'Color',[255,153,0]/256); end
if any(R_Aug>0); plot(ax4,nan,nan,'-','DisplayName','Auger','LineWidth',2,'Color',[255,25,80]/256); end
plot(ax4,nan,nan,'-','DisplayName','generation rate','LineWidth',2,'Color',[0,120,12]/256)
if any(Rl_SRH>0); plot(ax4,nan,nan,'xb','DisplayName','left surface SRH','LineWidth',1.5,'MarkerSize',10); end
if any(Rl_bim>0); plot(ax4,nan,nan,'ob','DisplayName','left surface bim','LineWidth',1.5,'MarkerSize',10); end
if any(Rr_SRH>0); plot(ax4,nan,nan,'xr','DisplayName','right surface SRH','LineWidth',1.5,'MarkerSize',10); end
if any(Rr_bim>0); plot(ax4,nan,nan,'or','DisplayName','right surface bim','LineWidth',1.5,'MarkerSize',10); end

% Add axis labels
xlabel(ax1,'applied voltage (V)');
ylabel(ax1,'current density (mAcm$^{-2}$)');
xlabel(ax2,'$x$ (nm)');
ylabel(ax2,'number density (m$^{-3}$)');
xlabel(ax3,'$x$ (nm)');
ylabel(ax3,'electric potential (V)');
xlabel(ax4,'$x$ (nm)');
ylabel(ax4,'recombination rate (m$^{-3}$s$^{-1}$)');

% Add timestamp
if log10(frame_times(end)-frame_times(1))<-1 % if time span is small use more precise timestamp
    prec = ceil(-log10(frame_times(end)-frame_times(1)))+3; % retain three significant figures
    str = ['$t =$ ' num2str(frame_times(1), prec) 's'];
else
    str = ['$t =$ ' datestr(seconds(frame_times(1)), 'HH:MM:SS.FFF')];
end
title(T,str,'FontSize',19,'Interpreter','latex');

% Add legends, position legend in centre of perovskite layer
leg2 = legend(ax2,'NumColumns',2,'Location','south');
leg2.BoxFace.ColorType = 'truecoloralpha';
leg2.BoxFace.ColorData = uint8(255*[1 1 1 0.5]'); % set transparancy
f = (-xE(1)+x(end)/2)/(-xE(1)+xH(end));
leg2.Position(1) = ax2.Position(1) + f*ax2.Position(3)-leg2.Position(3)/2;
leg4 = legend(ax4,'NumColumns',2,'Location','south');
leg4.BoxFace.ColorType = 'truecoloralpha';
leg4.BoxFace.ColorData = uint8(255*[1 1 1 0.5]'); % set transparancy

% Set axis properties
set(T, 'Padding', 'compact', 'TileSpacing', 'compact')
set(ax1,'YLim',YLims(1,:),'XLim',[floor(min(V)/0.1)*0.1,...
    ceil(max(V)/0.1)*0.1],'GridAlpha',0.4);
set(ax2,'YScale','log','YLim',YLims(2,:),'XLim',[xE(1) xH(end)]);
set(ax3,'YLim',YLims(3,:),'XLim',[xE(1) xH(end)]);
set(ax4,'YScale','log','YLim',YLims(4,:),'XLim',[x(1) x(end)]);
text(ax1,0.5,0.5,'J-V curve','Units','normalized','FontSize',30,'HorizontalAlignment','center')
text(ax2,0.5,0.5,'carrier densities','Units','normalized','FontSize',30,'HorizontalAlignment','center')
text(ax3,0.5,0.5,'electric potential','Units','normalized','FontSize',30,'HorizontalAlignment','center')
text(ax4,0.5,0.5,'recombination rates','Units','normalized','FontSize',30,'HorizontalAlignment','center')

% Output and close frame
frame = getframe(fig);
close(fignum);

end
