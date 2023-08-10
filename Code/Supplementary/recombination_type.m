function [Jd, Fi, nid] = recombination_type(type,params)
% A function to compute the parameters for plot_analytic_impedance and
% plot_analytic_current.

% Unpack parameters
[VT, Vbi, b, q, ni, beta, p0, tp, n0, tn, vpE, vnH] = ...
    struct2array(params, {'VT','Vbi','b','q','ni','beta','p0','tp', ...
                          'n0','tn','vpE','vnH'});

% Compute recombination type-dependent parameters
if strcmp(type,'Rb')
    nid = 1;
    Fi = @(Vs) Vs.V1+Vs.V2+Vs.V3+Vs.V4; % [V]
    Jd = q*b*beta*ni^2*exp(Vbi/VT)/10;  % [mA cm-2]
elseif strcmp(type,'Rp')
    nid = 2;
    Fi = @(Vs) Vs.V3+Vs.V4; % [V]
    Jd = q*b*p0/tp/10;      % [mA cm-2]
elseif strcmp(type,'Rn')
    nid = 2;
    Fi = @(Vs) Vs.V1+Vs.V2; % [V]
    Jd = q*b*n0/tn/10;      % [A cm-2]
elseif strcmp(type,'Rl')
    nid = 1;
    Fi = @(Vs) Vs.V2+Vs.V3+Vs.V4; % [V]
    Jd = q*vpE*p0/10;             % [mA cm-2]
elseif strcmp(type,'Rr')
    nid = 1;
    Fi = @(Vs) Vs.V1+Vs.V2+Vs.V3; % [V]
    Jd = q*vnH*n0/10;             % [mA cm-2]
else
    error('Unrecognised choice of recombination.');
end

% Check output
if isinf(Jd)
    Jd = 0;
end

end
