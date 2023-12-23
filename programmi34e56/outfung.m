function stop = outfung(optimValues,state)
global history
stop = false;
switch state
    case 'init'
        hold on; % initialized as empty
    case 'iter'
        history.fval = [history.fval; optimValues.localsolution.Fval];
        history.x = [history.x; optimValues.localsolution.X];
end
end

