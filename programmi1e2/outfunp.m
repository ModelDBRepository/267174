function [stop,options,optchanged] = myfun(optimvalues,options,flag)
global history
stop = false;
optchanged = 'true';
 
   switch flag
       case 'init'
           hold on
       case 'iter'
           % Concatenate current point and objective function
           % value with history. x must be a row vector.
           history.fval = [history.fval; optimvalues.fval];
           history.x = [history.x; optimvalues.x];
           save prova_history history
       case 'done'
           hold off
       otherwise
   end
end

