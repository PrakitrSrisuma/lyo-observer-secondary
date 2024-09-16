function outputs = PDE_Observer_Exp(t,y,T,ip)

% Important parameters and constants
m = ip.m;
dz = ip.dz;
k = ip.k;
Tb = min(ip.Tb0 + ip.r*t,ip.Tbmax);
ks = @(x) ip.A*exp(-ip.Ea/(ip.R*x));
q1 = k/(ip.rho*ip.Cp);
q2 = ip.rhod*ip.dHs/(ip.rho*ip.Cp);
K = ip.K;

% States
T_est = y(1:m);
cs_est = y(m+1:2*m);

% ODE
dcdt_est = zeros(m,1);
dTdt_est = zeros(m,1);

% Observer
for i = 1:m
    dcdt_est(i) = -ks(T_est(i))*cs_est(i);
end

% Heat transfer
for i = 2:m-1
    dTdt_est(i) = (q1/dz^2)*(T_est(i-1) - 2*T_est(i) + T_est(i+1)) + q2*dcdt_est(i);
end
dTdt_est(1) = 2*((q1/dz^2)*(T_est(2)-T_est(1))) + q2*dcdt_est(1);
dTdt_est(m) = 2*((q1/dz^2)*(T_est(m-1)-T_est(m))) + q2*dcdt_est(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T_est(m)-Tb);

switch ip.obs
case 'FullTemp'
LT = ip.LT*ones(m,m);
Lc = ip.Lc*ones(m,m);
y = T;
dTdt_est = dTdt_est + LT*(T_est-y);
dcdt_est = dcdt_est  + Lc*(T_est-y);

case 'OneTemp'
LT = ip.LT*ones(m,1);
Lc = ip.Lc*ones(m,1);
dTdt_est = dTdt_est + LT*(T_est(end)-T);
dcdt_est = dcdt_est  + Lc*(T_est(end)-T);

end

% Output
outputs = [dTdt_est;dcdt_est];

end