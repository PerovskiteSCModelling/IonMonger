function save_end_state(sol, filename)
% takes the final state of a solution file and prepares a vector
% for use as the initial condition in other simulations

disp('Saving end state as input vector')

params = sol.params;
[Vbi,VT,N,NE,NH,N0,dE,dH,kE,kH] = struct2array(params, {...
    'Vbi','VT','N','NE','NH','N0','dE','dH','kE','kH'});

% extract distributions at the end of solution file
P = sol.dstrbns.P(end,:)/N0;
phi = sol.dstrbns.phi(end,:)/VT;
n = sol.dstrbns.n(end,:)/(dE*kE);
p = sol.dstrbns.p(end,:)/(dH*kH);
phiE = sol.dstrbns.phiE(end,:)/VT;
nE = sol.dstrbns.nE(end,:)/dE;
phiH = sol.dstrbns.phiH(end,:)/VT;
pH = sol.dstrbns.pH(end,:)/dH;

% create column vector
u0 = [P';
    phi';
    n';
    p';
    phiE(1:end-1)';
    nE';
    phiH(2:end)';
    pH'];

Vapp = Vbi - VT*(phiE(1) - phiH(end));

inp_vec.u0 = u0;
inp_vec.Vapp = Vapp;

save(filename, 'inp_vec')

disp(['Input vector saved as ' filename '.mat'])

end