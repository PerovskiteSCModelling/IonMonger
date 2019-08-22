function vectors = create_vectors(params)
% This function creates the set of vectors required by the solver. The
% input is a structure containing the necessary parameters.

% ETL:   xE as NE+1 tanh-spaced points [-bE, 0   )
% Perov: x  as N +1 tanh-spaced points [  0, 1   ]
% HTL:   xH as NH+1 tanh-spaced points (  1, 1+bH]
% dxE, dx, dxH: as mesh spacings of xE, x, xH such that dx(i)=x(i+1)-x(i)
% Remove xE(end) and xH(1) as these points are x(1) and x(end) respectively

% Parameter input
[N, NE, NH, wE, wH] ...
    = struct2array(params,{'N','NE','NH','wE','wH'});

% Define spatial grids based on a "tanh" grid spacing
x = linspace(0,1,N+1)';
xE = linspace(-wE,0,NE+1)';
xH = linspace(1,1+wH,NH+1)';
st = params.st;
x = (tanh(st*(2*x-1))/tanh(st)+1)/2;
xE = wE*(tanh(st*(2*xE/wE+1))/tanh(st)-1)/2;
xH = 1+wH*(tanh(st*(2*(xH-1)/wH-1))/tanh(st)+1)/2;

% Define difference vectors
dx = diff(x);
dxE = diff(xE); % length NE
dxH = diff(xH); % length NH

% Output minimum and maximum spacing between grid points for each layer
if params.Verbose
    disp(['min dx  is ' num2str(min(dx )) ', max dx  is ' num2str(max(dx ))]);
    disp(['min dxE is ' num2str(min(dxE)) ', max dxE is ' num2str(max(dxE))]);
    disp(['min dxH is ' num2str(min(dxH)) ', max dxH is ' num2str(max(dxH))]);
end

% Redefine the spatial grids for the TLs to remove the points at the
% interfaces (which are included in the grid for the perovskite layer x)
xE = xE(1:end-1); % last point is first point of perovskite layer
xH = xH(2:end); % first point is last point of perovskite layer

% Package up vectors into a structure
vectors = struct('x',x,'dx',dx,'xE',xE,'dxE',dxE,'xH',xH,'dxH',dxH);

end
