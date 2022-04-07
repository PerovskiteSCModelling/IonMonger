function sol = compute_QFLs(sol)
% This script computes the quasi-Fermi levels and, if called as part of
% completion_tasks.m, these variables are added to the solution structure.

% Load parameters and variables
[VT, Ec, Ev, EcE, EvH, gc, gv, gcE, gvH] = ...
    struct2array(sol.params,{'VT','Ec','Ev','EcE','EvH','gc','gv', ...
                             'gcE','gvH'});
[phi, n, p, phiE, nE, phiH, pH] = ...
    struct2array(sol.dstrbns,{'phi','n','p','phiE','nE','phiH','pH'});

% Compute the quasi-Fermi levels and add them to the solution structure
sol.dstrbns.Efn  = -phi +Ec -VT*log(gc./n);
sol.dstrbns.Efp  = -phi +Ev +VT*log(gv./p);
sol.dstrbns.EfnE = -phiE+EcE-VT*log(gcE./nE);
sol.dstrbns.EfpH = -phiH+EvH+VT*log(gvH./pH);

end
