function outputs = get_input_data_exp2

% Secondary drying
input.mode = 'CFD';  % choose between 'HFD', 'CFD', 'MFD'
input.obs = 'FullTemp';  % choose between 'FullTemp' or 'OneTemp'
input.rho = 215;
input.rhod = 212.21;
input.Cp = 2590;  % [J/kgK]
input.Cpg = 1617;  % [J/kgK]
input.k = 0.217;
input.dHs = 2.68e6;  % heat of desoprtion [J/kg] 
input.H = 0.02;
input.d = 0.01;
input.h = 30;
input.r = .2/60;
input.Ea = 5700;
input.A = 1e-3;
input.R = 8.314;

% Microwave
input.Q = 85;
input.p1 = 3.73e-04;
input.p2 = 8.62e-03;
input.Qv_max = 6e3;
input.Qv_min = 1e2;
input.t_max = 6;
input.t_min = 6.5;

% Initial and boundary condition
input.T0 = 241.15;
input.Tb0 = 253.15; 
input.Tbmax = 313.15;
input.pin0 = 4;
input.pw0 = 1.07;
input.cs0 = .6415;
input.cs0_min = 0.0314;
input.cs0_max = 0.6415;
input.cfin = .01;

% State observer
input.LT = -1e-6;
input.Lc = 5e-7;

% Noise
input.noisemean = 0;  % in degree C
input.noisesd = 0;  % in degree C

% Controller gain
input.K = 1000;

% Numerics and Discretization
input.endtime = 12;  % final time (h) 
input.tol = 1e-6;  % tolerance
input.m = 20;  % number of nodes
input.dt = 60;  % (s)

% Export
outputs = input;

return