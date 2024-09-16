function [objective, termination, direction] = event_sublimation_completes(t,T,input)

objective = T(1)-1;  % stop when reaching sublimation temperature
termination = 1;  % terminate ode solvers 
direction = 0;  % both directions

end