function outputs = PDE_Observer_Scheduling(t,y,T_measured,tau,ip)

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

% Heat transfer
for i = 2:m-1
    dTdt(i) = (q1/dz^2)*(T(i-1) - 2*T(i) + T(i+1)) + q2*dcdt(i) + Qv ;
end
dTdt(1) = 2*((q1/dz^2)*(T(2)-T(1))) + q2*dcdt(1) + Qv ;
dTdt(m) = 2*((q1/dz^2)*(T(m-1)-T(m))) + q2*dcdt(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T(m)-Tb)+ Qv;

% Observer
for i = 1:m
    dcdt_est(i) = -ks(T_est(i))*cs_est(i);
end

% Heat transfer
for i = 2:m-1
    dTdt_est(i) = (q1/dz^2)*(T_est(i-1) - 2*T_est(i) + T_est(i+1)) + q2*dcdt_est(i)  + Qv;
end
dTdt_est(1) = 2*((q1/dz^2)*(T_est(2)-T_est(1))) + q2*dcdt_est(1) + Qv;
dTdt_est(m) = 2*((q1/dz^2)*(T_est(m-1)-T_est(m))) + q2*dcdt_est(m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*(T_est(m)-Tb) + Qv; 

if t >= tau*3600
    LT = ip.LT;
    Lc = 1e-7;
else
    LT = ip.LT;
    Lc = 5e-7;
end

switch ip.obs
case 'FullTemp'
LT = LT*ones(m,m);
Lc = Lc*ones(m,m);
y = zeros(m,1);
for i = 1:m
    y(i) = interp1(T_measured(:,end),T_measured(:,i),t);
end
dTdt_est = dTdt_est + LT*(T_est-y);
dcdt_est = dcdt_est  + Lc*(T_est-y);

case 'OneTemp'
LT = LT*ones(m,1);
Lc = Lc*ones(m,1);
dTdt_est = dTdt_est + LT*(mean(T_est)-mean(T));
dcdt_est = dcdt_est  + Lc*(mean(T_est)-mean(T));
% dTdt_est = dTdt_est + LT*(T_est(end)-T(end));
% dcdt_est = dcdt_est  + Lc*(T_est(end)-T(end));

end

% Output
outputs = [dTdt;dcdt;dTdt_est;dcdt_est];

end