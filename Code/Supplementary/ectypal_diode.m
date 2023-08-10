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
[time, tstar2t, jay, light, Vbi, psi2Vap, psi, VT, b, epsp] = ...
    struct2array(params, {'time','tstar2t','jay','light','Vbi','psi2Vap', ...
                          'psi','VT','b','epsp'});

% Compute dimensional voltage and then the time
V = psi2Vap(psi(time));
time = tstar2t(time);

% Compute the evolution of the Debye layer charge density and potentials
[Q, Vs] = ionic_charge(params);

% Compute the bulk electric field and its temporal derivative
Ebulk = (Vbi-V-Vs.V1-Vs.V2-Vs.V3-Vs.V4)/b;
dEdt = diff(Ebulk)./diff(time);
dEdt = ([dEdt(1),dEdt]+[dEdt,dEdt(end)])/2;

% Define a perovskite layer dimension
x = linspace(0,1,101)';

% Compute current loss for each bulk recombination mechanism
Rbulk = {'Rb','Rp','Rn'};
% Bimolecular
Jb = NaN(length(Rbulk),length(time));
[jr, Fi, nid] = recombination_type(Rbulk{1},params);
Jb(1,:) = jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);
% Hole-dominated SRH
[jr, Fi, nid] = recombination_type(Rbulk{2},params);
Jb(2,:) = jr.*exp(-Fi(Vs)/VT).*trapz(x,exp(-b*(1-x).*Ebulk/VT));
% Electron-dominated SRH
[jr, Fi, nid] = recombination_type(Rbulk{3},params);
Jb(3,:) = jr.*exp(-Fi(Vs)/VT).*trapz(x,exp(-b*x.*Ebulk/VT));

% Select dominant bulk mechanism
Jb(isinf(Jb)) = 0;
Jb = max(Jb);

% Compute current loss for ETL/perovskite interfacial recombination
[jr, Fi, nid] = recombination_type('Rl',params);
Jl = -jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);

% Compute current loss for ETL/perovskite interfacial recombination
[jr, Fi, nid] = recombination_type('Rr',params);
Jr = -jr.*exp(-(Fi(Vs)+b*Ebulk/nid)/VT);
    
% Compute the current densities
Jrec = Jb-Jl-Jr;        % total recombination current density [mA cm-2]
Jd = epsp*dEdt/10;           % displacement current density [mA cm-2]
J = jay*light(time)-Jrec+Jd; % total current density [mA cm-2]

% Package up solution
solution = struct('params',params, 'time',time, 'V',V', ...'Vres',Vres, ...
        'J',J', 'Jl',Jl', 'Jr',Jr', 'Jd',Jd', 'Q',Q', 'Vs',Vs);

% Save any impedance protocol parameters
if isfield(params,'impedance_protocol')
    solution.impedance_protocol = params.impedance_protocol;
end

end
