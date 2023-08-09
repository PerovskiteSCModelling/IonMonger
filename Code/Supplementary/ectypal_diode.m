function solution = ectypal_diode(params)
% A function to compute the current density from the ectypal diode
% equation for a given parameters structure. To produce a set of impedance
% simulations, use:
% for i = 1:length(sol)
%     params = sol(i).params;
%     params.impedance_protocol = sol(i).impedance_protocol;
%     sol_ec(i) = ectypal_diode(params);
% end

% Unpack parameters
[time, tstar2t, jay, light, Vbi, psi2Vap, psi, VT, b, epsp, q, N0] = ...
    struct2array(params, {'time','tstar2t','jay','light','Vbi','psi2Vap', ...
                          'psi','VT','b','epsp','q','N0'});

% Compute dimensional voltage and then the time
V = psi2Vap(psi(time));
time = tstar2t(time);

% Compute the evolution of the Debye layer charge density and potentials
[Q, Vs] = charge_evolution(params);
Q = sqrt(q*N0*epsp*VT)*Q; % [C m-2]

% Compute the bulk electric field and its temporal derivative
Ebulk = (Vbi-V-Vs.V1-Vs.V2-Vs.V3-Vs.V4)/b;
dEdt = diff(Ebulk)./diff(time);
dEdt = ([dEdt(1),dEdt]+[dEdt,dEdt(end)])/2;

% Compute current loss for each bulk recombination mechanism
Rbulk = {'Rb','Rp','Rn'};
Jb = NaN(length(Rbulk),length(time));
for i = 1:length(Rbulk)
    [jr, Fi, nid] = recombination_type(Rbulk{1},params);
    Jb(i,:) = jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);
end

% Compute current loss for ETL/perovskite interfacial recombination
[jr, Fi, nid] = recombination_type('Rl',params);
Jl = -jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);

% Compute current loss for ETL/perovskite interfacial recombination
[jr, Fi, nid] = recombination_type('Rr',params);
Jr = -jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);
    
% Sum up the recombination current densities
Jrec = max(Jb)-Jl-Jr;

% Compute the displacement current density
Jd = epsp*dEdt/10; % [mA cm-2]

% Compute the total current density
J = jay*light(time)-Jrec+Jd; % [mA cm-2]

% Package up solution
solution = struct('params',params, 'time',time, 'V',V', ...'Vres',Vres, ...
        'J',J', 'Jl',Jl', 'Jr',Jr', 'Jd',Jd', 'Q',Q', 'Vs',Vs);

% Save any impedance protocol parameters
if isfield(params,'impedance_protocol')
    solution.impedance_protocol = params.impedance_protocol;
end

end
