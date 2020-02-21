function dudt = RHS(t,u,psi,params,vectors,matrices,flag)
% This is the ODEFUN function for ode15s, containing the discretised
% right-hand sides of each of the DAEs. The additional inputs after the
% time and the vector of solution variables are the non-dimensional
% applied voltage and structures containing the necessary parameters,
% vectors and matrices for the calculation. The last input, flag, is an
% optional argument which is used to adapt the function, e.g. to simulate
% open-circuit conditions.

% Input parameters and arrays
[chi, delta, G, R, lambda, lam2, Rr, Rl, N, Kn, Kp, NE, lamE2, KE, kE, ...
    rE, NH, lamH2, KH, kH, rH, DI] ...
    = struct2array(params,{'chi','delta','G','R','lambda','lam2','Rr', ...
                           'Rl','N','Kn','Kp','NE','lamE2','KE','kE', ...
                           'rE','NH','lamH2','KH','kH','rH','DI'});
[x, dx, dxE, dxH] ...
    = struct2array(vectors,{'x','dx','dxE','dxH'});
[dudt, Av, AvE, AvH, Lo, LoE, LoH, Dx, DxE, DxH, NN, ddE, ddH] ...
    = struct2array(matrices,{'dudt','Av','AvE','AvH','Lo','LoE','LoH', ...
                             'Dx','DxE','DxH','NN','ddE','ddH'});

% Adjust for vectorisation
dudt = repmat(dudt,[1,size(u,2)]);

% Assign variable names
P   = u(1:N+1,:);
phi = u(N+2:2*N+2,:);
n   = u(2*N+3:3*N+3,:);
p   = u(3*N+4:4*N+4,:);
phiE = [u(4*N+5:4*N+NE+4,:); phi(1,:)];
nE   = [u(4*N+NE+5:4*N+2*NE+4,:); n(1,:)/kE];
phiH = [phi(end,:); u(4*N+2*NE+5:4*N+2*NE+NH+4,:)];
pH   = [p(end,:)/kH; u(4*N+2*NE+NH+5:4*N+2*NE+2*NH+4,:)];

% Compute variables (at the half points)
mE = Dx*phi; % negative electric field
mEE = DxE*phiE; % negative electric field in ETL
mEH = DxH*phiH; % negative electric field in HTL
FP = nnz(DI)*lambda*(Dx*P+mE.*(Av*P)); % negative anion vacancy flux
cd = NN-Lo*P+delta*(Lo*n-chi*Lo*p); % charge density
cdE = LoE*nE-ddE; % charge density in ETL
cdH = ddH-LoH*pH; % charge density in HTL
fn = Kn*(Dx*n-mE.*(Av*n)); % electron current
fnE = KE*(DxE*nE-mEE.*(AvE*nE)); % electron current in ETL
fp = Kp*(Dx*p+mE.*(Av*p)); % negative hole current
fpH = KH*(DxH*pH+mEH.*(AvH*pH)); % negative hole current in HTL
GR = G(Av*x,t)-R(Av*n,Av*p,Av*P); % generation-recombination

% P equation
dudt(1,:) = FP(1,:);
dudt(2:N,:) = FP(2:N,:)-FP(1:N-1,:);
dudt(N+1,:) = -FP(N,:);

% phi equation
dudt(N+2,:) = mE(1,:)-rE*mEE(end,:) ...
    -dx(1)*(1/2-P(1,:)/3-P(2,:)/6+delta*(n(1,:)/3+n(2,:)/6-chi*(p(1,:)/3+p(2,:)/6)))/lam2 ...
    -rE*dxE(end)*(nE(end-1,:)/6+nE(end,:)/3-1/2)/lamE2; % continuity
dudt(N+3:2*N+1,:) = mE(2:N,:)-mE(1:N-1,:)-cd/lam2;
dudt(2*N+2,:) = rH*mEH(1,:)-mE(end,:) ...
    -dx(end)*(1/2-P(end-1,:)/6-P(end,:)/3+delta*(n(end-1,:)/6+n(end,:)/3-chi*(p(end-1,:)/6+p(end,:)/3)))/lam2 ...
    -rH*dxH(1)*(1/2-pH(1,:)/3-pH(2,:)/6)/lamH2; % continuity

% n equation
dudt(2*N+3,:) = fn(1,:)-fnE(end,:)+(dx(1).*GR(1,:))/2-Rl(nE(end,:),p(1,:)); % continuity
dudt(2*N+4:3*N+2,:) = fn(2:N,:)-fn(1:N-1,:)+(dx(2:N).*GR(2:N,:)+dx(1:N-1).*GR(1:N-1,:))/2;
dudt(3*N+3,:) = -fn(N,:)+dx(N)*GR(end,:)/2-Rr(n(N+1,:),pH(1,:));

% p equation
dudt(3*N+4,:) = fp(1,:)+dx(1)*GR(1,:)/2-Rl(nE(end,:),p(1,:));
dudt(3*N+5:4*N+3,:) = fp(2:N,:)-fp(1:N-1,:)+(dx(2:N).*GR(2:N,:)+dx(1:N-1).*GR(1:N-1,:))/2;
dudt(4*N+4,:) = fpH(1,:)-fp(end,:)+(dx(end)*GR(end,:))/2-Rr(n(N+1,:),pH(1,:)); % continuity

% phiE equation
dudt(4*N+5,:) = phiE(1,:)-psi(t);
dudt(4*N+6:4*N+NE+4,:) = mEE(2:NE,:)-mEE(1:NE-1,:)-cdE/lamE2;

% nE equation
dudt(4*N+NE+5,:) = nE(1,:)-1;
dudt(4*N+NE+6:4*N+2*NE+4,:) = fnE(2:NE,:)-fnE(1:NE-1,:);

% phiH equation
dudt(4*N+2*NE+5:4*N+2*NE+NH+3,:) = mEH(2:NH,:)-mEH(1:NH-1,:)-cdH/lamH2;
dudt(4*N+2*NE+NH+4,:) = phiH(end,:)+psi(t);

% pH equation
dudt(4*N+2*NE+NH+5:4*N+2*NE+2*NH+3,:) = fpH(2:NH,:)-fpH(1:NH-1,:);
dudt(4*N+2*NE+2*NH+4,:) = pH(end,:)-1;

% Perform any additional step requested by the optional input argument flag
if nargin>6
    if strcmp(flag,'none')
        % Do nothing else
    elseif strcmp(flag,'open-circuit')
        % Overwrite the entries for the potential at each contact
        % Zero current boundary condition:
        dudt(4*N+5,:) = fnE(1,:);
        % Symmetric values of the potential at the contacts
        dudt(4*N+2*NE+NH+4,:) = phiE(1,:)+phiH(end,:);
    elseif strcmp(flag,'init')
        % Overwrite right-hand BC to ensure conservation of ion vacancies
        dudt(N+1,:) = trapz(x,P)-1;
    elseif strcmp(flag,'findVoc')
        % Overwrite the entries for the potential at each contact
        % Zero current boundary condition:
        dudt(4*N+5,:) = fnE(1,:);
        % Symmetric values of the potential at the contacts
        dudt(4*N+2*NE+NH+4,:) = phiE(1,:)+phiH(end,:);
        % Overwrite right-hand BC to ensure conservation of ion vacancies
        dudt(N+1,:) = trapz(x,P)-1;
    else
        error('The optional input argument is not recognised.');
    end
end

end
