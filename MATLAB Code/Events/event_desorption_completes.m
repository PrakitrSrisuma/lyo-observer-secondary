function [objective, termination, direction] = event_desorption_completes(t,T,input)

objective = T(end)-input.cfin;  % stop when reaching sublimation temperature
termination = 1;  % terminate ode solvers 
direction = 0;  % both directions

end