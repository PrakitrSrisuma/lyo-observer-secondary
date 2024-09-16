function outputs = cal_timeconstant(ip)

m = ip.m;
dx = 1e-3;
T0 = 0.5*(ip.T0+ip.Tbmax)*ones(m,1);
cs0 = 0.5*(ip.cs0)*ones(m,1);

y0 = [T0;cs0];
A = PDE_Jacobian_Analytical(y0,ip);
A2 = PDE_Jacobian_Numerical(y0,ip,dx);

LT = ip.LT*ones(m,m);
Lc = ip.Lc*ones(m,m);
L = [LT;Lc];
C = [eye(m),zeros(m)];


B = A+L*C;
B2 = A2+L*C;
[V,D] = eig(B);
U = inv(V);

e = eig(B);
e2 = eig(B2);
e0 = eig(A);
tau = abs(real(e(m+1)));
tau = (1/tau)/3600;
tau2 = abs(real(e2(m+1)));
tau2 = (1/tau2)/3600;

outputs.tau = tau;
outputs.tau_0 = 4*tau;
outputs.tau2 = tau2;
outputs.tau2_0 = 4*tau2;
outputs.U = U;
outputs.V = V;
outputs.y0 = y0;
outputs.eig = e;

return