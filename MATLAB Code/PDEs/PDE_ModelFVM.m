function outputs = PDE_ModelFVM(t,y,ip)

% Important parameters and constants
m = ip.m;
dz = ip.dz;
k = ip.k;
Tb = min(ip.Tb0 + ip.r*t,ip.Tbmax);
ks = @(x) ip.A*exp(-ip.Ea/(ip.R*x));
q1 = k/(ip.rho*ip.Cp);
q2 = ip.rhod*ip.dHs/(ip.rho*ip.Cp);
Qv = interp1(ip.Qv(:,1),ip.Qv(:,2),t);
Qv = Qv/(ip.rho*ip.Cp);

% States
T = y(1:m);
cs = y(m+1:2*m);

% ODE
dcdt = zeros(m,1);
dTdt = zeros(m,1);

% Desorption
for i = 1:m
    dcdt(i) = -ks(T(i))*cs(i);
    %dcdt(i) = -ip.ks0*cs(i);
end

ks_check = ks(T(i));
% Heat transfer
for i = 2:m-1
    dTdt(i) = (q1/dz^2)*(T(i-1) - 2*T(i) + T(i+1)) + q2*dcdt(i) + Qv ;
end
dTdt(1) = 2*((q1/dz^2)*(T(2)-T(1))) + q2*dcdt(1) + Qv ;
dTdt(m) = 2*((q1/dz^2)*(T(m-1)-T(m))) + q2*dcdt(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T(m)-Tb)+ Qv;

% Output
outputs = [dTdt;dcdt];

end