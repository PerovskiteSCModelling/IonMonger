function matrices = create_matrices(params,vectors)
% This function creates the set of matrices required by the solver. The
% inputs are structures containing the necessary parameters and vectors.

% Parameter input
[N, NE, NH] = struct2array(params,{'N','NE','NH'});
[dx, dxE, dxH] = struct2array(vectors,{'dx','dxE','dxH'});

% Preallocation for dudt
dudt = spalloc(4*N+2*NE+2*NH+4,1,4*N+2*NE+2*NH+4);

% Define averaging matrices
Av = gallery('tridiag',N+1,0,1,1)/2; Av = Av(1:N,1:N+1);
AvE = gallery('tridiag',NE+1,0,1,1)/2; AvE = AvE(1:NE,1:NE+1);
AvH = gallery('tridiag',NH+1,0,1,1)/2; AvH = AvH(1:NH,1:NH+1);
Lo = gallery('tridiag',N+1,[dx(1:end-1)/6; 0],[0; (dx(1:end-1)+dx(2:end))/3; 0],[0; dx(2:end)/6]);
Lo = Lo(2:N,1:N+1);
LoE = gallery('tridiag',NE+1,[dxE(1:end-1)/6; 0],[0; (dxE(1:end-1)+dxE(2:end))/3; 0],[0; dxE(2:end)/6]);
LoE = LoE(2:NE,1:NE+1);
LoH = gallery('tridiag',NH+1,[dxH(1:end-1)/6; 0],[0; (dxH(1:end-1)+dxH(2:end))/3; 0],[0; dxH(2:end)/6]);
LoH = LoH(2:NH,1:NH+1);

% Define differencing matrices
Dx = gallery('tridiag',N+1,0./dx,[-1./dx; 0],1./dx); Dx = Dx(1:N,1:N+1);
DxE = gallery('tridiag',NE+1,0./dxE,[-1./dxE; 0],1./dxE); DxE = DxE(1:NE,1:NE+1);
DxH = gallery('tridiag',NH+1,0./dxH,[-1./dxH; 0],1./dxH); DxH = DxH(1:NH,1:NH+1);

% Define vectors for constant cation vacancy and doping densities
NN = (dx(2:end)+dx(1:end-1))/2;
ddE = (dxE(2:end)+dxE(1:end-1))/2;
ddH = (dxH(2:end)+dxH(1:end-1))/2;

% Package up matrices into a structure
matrices = struct('dudt',dudt,'Av',Av,'AvE',AvE,'AvH',AvH, ...
                              'Lo',Lo,'LoE',LoE,'LoH',LoH, ...
                              'Dx',Dx,'DxE',DxE,'DxH',DxH, ...
                              'NN',NN,'ddE',ddE,'ddH',ddH);

end
