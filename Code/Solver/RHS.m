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
    rE, NH, lamH2, KH, kH, rH, DI, nc, pc, ARs, Rsp, pbi, nonlinear, lim, SEinv, ...
    omegE, SHinv, omegH, AE, AH, phidisp] ...
    = struct2array(params,{'chi','delta','G','R','lambda','lam2','Rr', ...
                           'Rl','N','Kn','Kp','NE','lamE2','KE','kE', ...
                           'rE','NH','lamH2','KH','kH','rH','DI','nc', ...
                           'pc','ARs','Rsp','pbi', 'nonlinear', 'lim','SEinv', ...
                           'omegE', 'SHinv', 'omegH', 'AE', 'AH', 'phidisp'});
[x, dx, dxE, dxH] = struct2array(vectors,{'x','dx','dxE','dxH'});
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
nE   = u(4*N+NE+5:4*N+2*NE+5,:);
phiH = [phi(end,:); u(4*N+2*NE+6:4*N+2*NE+NH+5,:)];
pH   = u(4*N+2*NE+NH+6:4*N+2*NE+2*NH+6,:);

% Compute variables (at the half points)
mE = Dx*phi; % negative electric field
mEE = DxE*phiE; % negative electric field in ETL
mEH = DxH*phiH; % negative electric field in HTL

if any(lim) && strcmp(nonlinear,'Drift')
    PAP = lim*(2*Av*(P.^2)+P(1:N,:).*P(2:N+1,:))/3;
    PD = Dx*P;
elseif any(lim) && strcmp(nonlinear,'Diffusion')
    PAP = 0;
    PD = -1/lim*Dx*log(1-lim*P);
else
  PAP = 0;
	PD = Dx*P;
end
FP = nnz(DI)*lambda*(PD+mE.*(Av*P-PAP)); % negative anion vacancy flux
cd = NN-Lo*P+delta*(Lo*n-chi*Lo*p); % charge density
cdE = LoE*nE-ddE; % charge density in ETL
cdH = ddH-LoH*pH; % charge density in HTL
fn = Kn*(Dx*n-mE.*(Av*n)); % electron current
fnE = kE*KE*(AvE*nE).*(DxE*(SEinv(omegE*nE)-phiE)); % electron current in ETL
fp = Kp*(Dx*p+mE.*(Av*p)); % negative hole current
fpH = kH*KH*(AvH*pH).*(DxH*(SHinv(omegH*pH)+phiH)); % negative hole current in HTL
GR = G(Av*x,t)-R(Av*n,Av*p,Av*P); % generation-recombination

% P equation
dudt(1,:) = FP(1,:);
dudt(2:N,:) = FP(2:N,:)-FP(1:N-1,:);
dudt(N+1,:) = -FP(N,:);

% phi equation
dudt(N+2,:) = mE(1,:)-rE*mEE(end,:) ...
    -dx(1)*(1/2-P(1,:)/3-P(2,:)/6+delta*(n(1,:)/3+n(2,:)/6-chi*(p(1,:)/3+...
    p(2,:)/6)))/lam2-...
    rE*dxE(end)*(nE(end-1,:)/6+nE(end,:)/3-1/2)/lamE2; % continuity
dudt(N+3:2*N+1,:) = mE(2:N,:)-mE(1:N-1,:)-cd/lam2;
dudt(2*N+2,:) = rH*mEH(1,:)-mE(end,:) ...
    -dx(end)*(1/2-P(end-1,:)/6-P(end,:)/3+delta*(n(end-1,:)/6+n(end,:)/3-...
    chi*(p(end-1,:)/6+p(end,:)/3)))/lam2 ...
    -rH*dxH(1)*(1/2-pH(1,:)/3-pH(2,:)/6)/lamH2; % continuity

% n equation
dudt(2*N+3,:) = kE*fn(1,:)-fnE(end,:)+kE*(dx(1).*GR(1,:)/2-Rl(n(1,:),p(1,:))); % continuity
dudt(2*N+4:3*N+2,:) = fn(2:N,:)-fn(1:N-1,:)+(dx(2:N).*GR(2:N,:)+dx(1:N-1).*GR(1:N-1,:))/2;
dudt(3*N+3,:) = -fn(N,:)+dx(N)*GR(end,:)/2-Rr(n(N+1,:),p(N+1,:));

% p equation
dudt(3*N+4,:) = fp(1,:)+dx(1)*GR(1,:)/2-Rl(n(1,:),p(1,:));
dudt(3*N+5:4*N+3,:) = fp(2:N,:)-fp(1:N-1,:)+(dx(2:N).*GR(2:N,:)+dx(1:N-1).*GR(1:N-1,:))/2;
dudt(4*N+4,:) = fpH(1,:)-kH*fp(end,:)+kH*(dx(end)*GR(end,:)/2-Rr(n(N+1,:),p(N+1,:))); % continuity

% phiE equation
dudt(4*N+5,:) = phiE(1,:)-psi(t)-phidisp;
dudt(4*N+6:4*N+NE+4,:) = mEE(2:NE,:)-mEE(1:NE-1,:)-cdE/lamE2;

% nE equation
dudt(4*N+NE+5,:) = nE(1,:)-nc;
dudt(4*N+NE+6:4*N+2*NE+4,:) = fnE(2:NE,:)-fnE(1:NE-1,:);
dudt(4*N+2*NE+5,:) = n(1,:) - exp(SEinv(omegE*nE(end,:))-SEinv(omegE));

% phiH equation
dudt(4*N+2*NE+6:4*N+2*NE+NH+4,:) = mEH(2:NH,:)-mEH(1:NH-1,:)-cdH/lamH2;
dudt(4*N+2*NE+NH+5,:) = phiH(end,:)+psi(t) ...
                        +(fpH(end,:)/kH*ARs+Rsp*(pbi-2*psi(t)))/(1+Rsp)-phidisp;
% Note that the last term in this BC accounts for any parasitic resistance
% (neglecting the displacement current, which should be small at the contacts)

% pH equation
dudt(4*N+2*NE+NH+6,:) = p(end,:) - exp(SHinv(omegH*pH(1,:))-SHinv(omegH));
dudt(4*N+2*NE+NH+7:4*N+2*NE+2*NH+5,:) = fpH(2:NH,:)-fpH(1:NH-1,:);
dudt(4*N+2*NE+2*NH+6,:) = pH(end,:)-pc;

% Perform any additional step requested by the optional input argument flag
if nargin>6
    if strcmp(flag,'none')
        % Do nothing else
    elseif strcmp(flag,'open-circuit')
        % Overwrite the entries for the potential at each contact
        % Zero current boundary condition:
        dudt(4*N+5,:) = fnE(1,:);
        % Symmetric values of the potential at the contacts
        dudt(4*N+2*NE+NH+5,:) = phiE(1,:)+phiH(end,:)-2*phidisp;
    elseif strcmp(flag,'init')
        % Overwrite right-hand BC to ensure conservation of ion vacancies
        dudt(N+1,:) = trapz(x,P)-1;
    elseif strcmp(flag,'findVoc')
        % Overwrite the entries for the potential at each contact
        % Zero current boundary condition:
        dudt(4*N+5,:) = fnE(1,:);
        % Symmetric values of the potential at the contacts
        dudt(4*N+2*NE+NH+5,:) = phiE(1,:)+phiH(end,:)-2*phidisp;
        % Overwrite right-hand BC to ensure conservation of ion vacancies
        dudt(N+1,:) = trapz(x,P)-1;
    else
        error('The optional input argument is not recognised.');
    end
end

end
