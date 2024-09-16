function outputs = PDE_Jacobian_Analytical(y,ip)

% Important parameters and constants
m = ip.m;
dz = ip.dz;
k = ip.k;
ks = @(x) ip.A*exp(-ip.Ea/(ip.R*x));
q1 = k/(ip.rho*ip.Cp);
q2 = ip.rhod*ip.dHs/(ip.rho*ip.Cp);

% States
T = y(1:m);
cs = y(m+1:2*m);

% Jacobian
J = zeros(2*m);

% Desorption
for i = 2:m-1
    J(i,i) = -2*q1/dz^2 - q2*cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2));
    J(i,i+m) = -q2*ks(T(i));
    J(i,i+1) = q1/dz^2;
    J(i,i-1) = q1/dz^2;

    J(i+m,i+m) = -ks(T(i));
    J(i+m,i) = -cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2));
end

for i = 1
    J(i,i) = -2*q1/dz^2 - q2*cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2));
    J(i,i+m) = -q2*ks(T(i));
    J(i,i+1) = 2*q1/dz^2;

    J(i+m,i+m) = -ks(T(i));
    J(i+m,i) = -cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2));   
end

for i = m
    J(i,i) = -2*q1/dz^2 - q2*cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2)) - 2*(ip.h/(ip.rho*ip.Cp*dz));
    J(i,i+m) = -q2*ks(T(i));
    J(i,i-1) = 2*q1/dz^2;

    J(i+m,i+m) = -ks(T(i));
    J(i+m,i) = -cs(i)*ks(T(i))*(ip.Ea/(ip.R*T(i)^2)); 
end

% Output
outputs = J;

end