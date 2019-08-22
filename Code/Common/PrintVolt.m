function status = PrintVolt(~,u,flag,params)
% An ode15s OutputFcn which can be used to print the applied voltage at each step

N = params.N;
if isempty(flag)
    disp(['V = ',num2str(double(params.psi2Vap(u(4*N+5)))), ' V']);
elseif isequal(flag,'init')
    disp(['Starting, V = ',num2str(double(params.psi2Vap(u(4*N+5)))),' V']);
else
    disp('Finished.'); fprintf('\n\n');
end
status = 0;

end