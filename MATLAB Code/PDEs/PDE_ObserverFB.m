function outputs = PDE_ObserverFB(t,y,ip)

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
T = y(1:m);
cs = y(m+1:2*m);
T_est = y(2*m+1:3*m);
cs_est = y(3*m+1:end);

% ODE
dcdt = zeros(m,1);
dTdt = zeros(m,1);
dcdt_est = zeros(m,1);
dTdt_est = zeros(m,1);

% Desorption
for i = 1:m
    dcdt(i) = -ks(T(i))*cs(i);
    %dcdt(i) = -ip.ks0*cs(i);
end

Tmax = max(T);
Tmax_est = max(T_est);

% Heat transfer
for i = 2:m-1
    dTdt(i) = (q1/dz^2)*(T(i-1) - 2*T(i) + T(i+1)) + q2*dcdt(i) + (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax);
end
dTdt(1) = 2*((q1/dz^2)*(T(2)-T(1))) + q2*dcdt(1) + (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax);
dTdt(m) = 2*((q1/dz^2)*(T(m-1)-T(m))) + q2*dcdt(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T(m)-Tb)+ (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax) ;

% Observer
for i = 1:m
    dcdt_est(i) = -ks(T_est(i))*cs_est(i);
end

% Heat transfer
for i = 2:m-1
    dTdt_est(i) = (q1/dz^2)*(T_est(i-1) - 2*T_est(i) + T_est(i+1)) + q2*dcdt_est(i)  + (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax_est);
end
dTdt_est(1) = 2*((q1/dz^2)*(T_est(2)-T_est(1))) + q2*dcdt_est(1) + (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax_est);
dTdt_est(m) = 2*((q1/dz^2)*(T_est(m-1)-T_est(m))) + q2*dcdt_est(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T_est(m)-Tb)+ (K/(ip.rho*ip.Cp))*(ip.Tbmax-Tmax_est) ;

switch ip.obs
case 'FullTemp'
LT = ip.LT*ones(m,m);
Lc = ip.Lc*ones(m,m);
y = T + interp1(ip.noise(:,2),ip.noise(:,1),t);
dTdt_est = dTdt_est + LT*(T_est-y);
dcdt_est = dcdt_est  + Lc*(T_est-y);

case 'OneTemp'
LT = ip.LT*ones(m,1);
Lc = ip.Lc*ones(m,1);
dTdt_est = dTdt_est + LT*(mean(T_est)-mean(T));
dcdt_est = dcdt_est  + Lc*(mean(T_est)-mean(T));

end

% Output
outputs = [dTdt;dcdt;dTdt_est;dcdt_est];

end