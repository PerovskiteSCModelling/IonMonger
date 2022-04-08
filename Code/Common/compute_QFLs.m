function sol = compute_QFLs(sol)
% This script computes the quasi-Fermi levels and, if called as part of
% completion_tasks.m, these variables are added to the solution structure.

% Load parameters and variables
[VT, Ec, Ev, gc, gv, nE2EfE, pH2EfH]  = ...
    struct2array(sol.params,{'VT','Ec','Ev','gc','gv','nE2EfE','pH2EfH'});
[EcE, EvH, gcE, gvH] = struct2array(sol.params,{'EcE','EvH','gcE','gvH'});
[phi, n, p, phiE, nE, phiH, pH] = ...
    struct2array(sol.dstrbns,{'phi','n','p','phiE','nE','phiH','pH'});

% Compute the quasi-Fermi levels and add them to the solution structure
sol.dstrbns.Efn  = -phi +Ec -VT*log(gc./n);
sol.dstrbns.Efp  = -phi +Ev +VT*log(gv./p);
if isfield(sol.params,'SE') % Check if solution is from v1 or v2
    sol.dstrbns.EfnE = -phiE+nE2EfE(nE); % use stats function
else
    sol.dstrbns.EfnE = -phiE+EcE-VT*log(gcE./nE); % assume Boltzmann
end
if isfield(sol.params,'SH') % Check if solution is from v1 or v2
    sol.dstrbns.EfpH = -phiH+pH2EfH(pH); % use stats function
else
    sol.dstrbns.EfpH = -phiH+EvH+VT*log(gvH./pH); % assume Boltzmann
end

end
