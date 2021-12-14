function status = PrintTime(t,~,flag,params)
% An ode15s OutputFcn which can be used to print the time at each step

if isempty(flag)
    disp(['t = ',num2str(double(params.Tion*t(end))), ' s']);
elseif isequal(flag,'init')
    disp(['Starting, t = ',num2str(double(params.Tion*t(1))),' s']);
else
    disp('Finished.'); fprintf('\n\n');
end
status = 0;

end