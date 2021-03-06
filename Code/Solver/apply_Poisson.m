function sol_init = apply_Poisson(sol_init,params,vectors,matrices)
% This function applies the algebraic equations corresponding to Poisson's
% equation for the electric potential to a vector of the solution variables
% (sol_init) in order to ensure that Poisson's equation is satisfied as
% precisely as can be achieved by mldivide. The other inputs are structures
% containing the necessary parameters, vectors and matrices for the
% computation.

% Parameter input
[delta, chi, lam2, lamE2, lamH2, N, NE, NH, kH, kE, rE, rH] ...
    = struct2array(params, {'delta','chi','lam2','lamE2','lamH2', ...
    'N','NE','NH','kH','kE','rE','rH'});
[dx, dxE, dxH] = struct2array(vectors, {'dx','dxE','dxH'});
[Lo, LoE, LoH, NN, ddE, ddH] = struct2array(matrices,{'Lo','LoE','LoH','NN','ddE','ddH'});

% Algebraic equations corresponding to Poisson's equation
A = gallery('tridiag',N+NE+NH-1, ...
    [1./dxE(2:end-1); rE./dxE(end); 1./dx; 1./dxH(1:end-1)], ...
    [-(1./dxE(1:end-1)+1./dxE(2:end)); -(rE./dxE(end)+1./dx(1)); ...
    -(1./dx(1:end-1)+1./dx(2:end)); -(1./dx(end)+rH./dxH(1)); ...
    -(1./dxH(1:end-1)+1./dxH(2:end))], ...
    [1./dxE(2:end); 1./dx; rH./dxH(1); 1./dxH(2:end-1)]);
B(1:NE-1,1) = (LoE*[sol_init(4*N+NE+5:4*N+2*NE+4); sol_init(2*N+3)/kE] ...
        -ddE)/lamE2; % charge density in ETL
B(1) = B(1)-sol_init(4*N+5)/dxE(1); % to take account of contact potential
B(NE,1) = dx(1)*(1/2-sol_init(1)/3-sol_init(2)/6 ...
        +delta*(sol_init(2*N+3)/3+sol_init(2*N+4)/6 ...
        -chi*(sol_init(3*N+4)/3+sol_init(3*N+5)/6)))/lam2 ...
        +rE*dxE(end)*(sol_init(4*N+2*NE+4)/6 ...
        +sol_init(2*N+3)/kE/3-1/2)/lamE2; % continuity
B(NE+1:N+NE-1,1) = (NN-Lo*sol_init(1:N+1)+delta*(Lo*sol_init(2*N+3:3*N+3) ...
        -chi*Lo*sol_init(3*N+4:4*N+4)))/lam2; % charge density
B(N+NE,1) = dx(end)*(1/2-sol_init(N)/6-sol_init(N+1)/3 ...
        +delta*(sol_init(3*N+2)/6+sol_init(3*N+3)/3 ...
        -chi*(sol_init(4*N+3)/6+sol_init(4*N+4)/3)))/lam2 ...
        +rH*dxH(1)*(1/2-sol_init(4*N+4)/kH/3 ...
        -sol_init(4*N+2*NE+NH+5)/6)/lamH2; % continuity
B(N+NE+1:N+NE+NH-1,1) = (ddH-LoH*[sol_init(4*N+4)/kH; ...
        sol_init(4*N+2*NE+NH+5:4*N+2*NE+2*NH+4)]) ...
        /lamH2; % charge density in HTL
B(N+NE+NH-1) = B(N+NE+NH-1)-sol_init(4*N+2*NE+NH+4)/dxH(end); % contact
sol_init([4*N+6:4*N+NE+4,N+2:2*N+2,4*N+2*NE+5:4*N+2*NE+NH+3]) = mldivide(A,B);

end