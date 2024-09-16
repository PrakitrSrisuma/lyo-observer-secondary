function outputs = PDE_Jacobian_Numerical(y,ip,dx)

% Important parameters and constants
m = ip.m;
dz = ip.dz;
k = ip.k;
ks = @(x) ip.A*exp(-ip.Ea/(ip.R*x));
q1 = k/(ip.rho*ip.Cp);
q2 = ip.rhod*ip.dHs/(ip.rho*ip.Cp);


J = zeros(2*m);

for i = 1
    J(i,i) = (1/dx)*(((2*q1/dz^2)*(y(i+1) - (y(i)+dx)) - q2*ks(y(i)+dx)*y(i+m))- ((2*q1/dz^2)*(y(i+1) - y(i)) - q2*ks(y(i))*y(i+m)));
    J(i,i+1) = (1/dx)*(((2*q1/dz^2)*(y(i+1)+dx - y(i)) - q2*ks(y(i))*y(i+m))- ((2*q1/dz^2)*(y(i+1) - y(i)) - q2*ks(y(i))*y(i+m)));
    J(i,i+m) = (1/dx)*(((2*q1/dz^2)*(y(i+1) - y(i)) - q2*ks(y(i))*(y(i+m)+dx))- ((2*q1/dz^2)*(y(i+1) - y(i)) - q2*ks(y(i))*y(i+m)));
end

for i = 2:m-1
    J(i,i-1) = (1/dx)*(((q1/dz^2)*(y(i-1)+dx - 2*y(i) + y(i+1)) - q2*ks(y(i))*y(i+m)) - ((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)) - q2*ks(y(i))*y(i+m)));
    J(i,i) = (1/dx)*(((q1/dz^2)*(y(i-1) - 2*(y(i)+dx) + y(i+1)) - q2*ks(y(i)+dx)*y(i+m))- ((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)) - q2*ks(y(i))*y(i+m)));
    J(i,i+1) = (1/dx)*(((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)+dx) - q2*ks(y(i))*y(i+m))-((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)) - q2*ks(y(i))*y(i+m)));
    J(i,i+m) = (1/dx)*(((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)) - q2*ks(y(i))*(y(i+m)+dx))- ((q1/dz^2)*(y(i-1) - 2*y(i) + y(i+1)) - q2*ks(y(i))*y(i+m)));
end

for i = m
    J(i,i-1) = (1/dx)*(((2*q1/dz^2)*(y(i-1)+dx - y(i)) - q2*ks(y(i))*y(i+m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*y(i))- ((2*q1/dz^2)*(y(i-1) - y(i)) - q2*ks(y(i))*y(i+m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*y(i)));
    J(i,i) = (1/dx)*(((2*q1/dz^2)*(y(i-1) - (y(i)+dx)) - q2*ks(y(i)+dx)*y(i+m)- 2*(ip.h/(ip.rho*ip.Cp*dz))*(y(i)+dx))- ((2*q1/dz^2)*(y(i-1) - y(i)) - q2*ks(y(i))*y(i+m)- 2*(ip.h/(ip.rho*ip.Cp*dz))*y(i)));
    J(i,i+m) = (1/dx)*(((2*q1/dz^2)*(y(i-1) - y(i)) - q2*ks(y(i))*(y(i+m)+dx) - 2*(ip.h/(ip.rho*ip.Cp*dz))*y(i))- ((2*q1/dz^2)*(y(i-1) - y(i)) - q2*ks(y(i))*y(i+m) - 2*(ip.h/(ip.rho*ip.Cp*dz))*y(i)));
end

for i = m+1:2*m
    J(i,i) = (1/dx)*(-ks(y(i-m))*(y(i)+dx) + ks(y(i-m))*(y(i)));
    J(i,i-m) = (1/dx)*(-ks(y(i-m)+dx)*(y(i)) + ks(y(i-m))*(y(i)));
end

% Output
outputs = J;

end