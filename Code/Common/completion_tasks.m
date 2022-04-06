function sol = completion_tasks(sol)
% This script is called by numericalsolver and can contain any additional
% tasks to be performed after the (dimensional) solution has been found.
% The only input is the solution structure.

% Compute the quasi-Fermi levels
% sol = compute_QFLs(sol);

if sol.params.Verbose
    
    % Calculate the absolute percentage error in the conservation of ions
    [N0, b] = struct2array(sol.params,{'N0','b'});
    x = struct2array(sol.vectors,{'x'});
    P = struct2array(sol.dstrbns,{'P'});
    fprintf('The absolute percentage error in the anion mass conservation is %0.2g%% \n', ...
        abs(trapz(x*1e-9,N0-P(end,:)))/(b*N0));
    
    % Output the steady-state Voc if computed
    if sol.params.findVoc
        phiE = sol.dstrbns.phiE/sol.params.VT;
        fprintf('The Voc was found to be %0.5g V \n', ...
            sol.params.psi2Vap(phiE(1,1)));
    end
    
end

end
