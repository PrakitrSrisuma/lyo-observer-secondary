function outputs = input_processing(input)

% Secondary drying
input.V = pi*input.d^2*input.H/4;  % volume of one vial (m3)
input.Dg0 = @(T) 0.000143931*(T^3*(1/input.Mw + 1/input.Min))^0.5;
input.Dg = @(Dg0,P) Dg0/P;
input.k1 = @(Dg0,Kw,Kg,P) input.C2*Dg0*Kw/(input.C2*Dg0+Kg*P);
input.k2 = @(Dg0,Kw,Kin,Kg,mug,P) Kin*Kw/(input.C2*Dg0+Kg*P) + input.C0/mug;
input.k3 = @(Dg0,Kin,Kg,P) input.C2*Dg0*Kin/(input.C2*Dg0+Kg*P);
input.Kw = @(T) input.C1/(input.R*T/input.Mw)^0.5;
input.Kin = @(T) input.C1/(input.R*T/input.Min)^0.5;
input.p0 = input.pw0 + input.pin0;

switch input.mode
case 'MFD'
input.h = 0;
input.tv = [0; input.t_max; input.t_min; input.endtime]*3600;
input.Qv = [input.Qv_max ; input.Qv_max ; input.Qv_min; input.Qv_min];
input.Qv = [input.tv,input.Qv];
case 'CFD'
input.tv = [0; input.endtime]*3600;
input.Qv = [0;0];
input.Qv = [input.tv,input.Qv];
input.K = 0;
case 'HFD'
input.tv = [0; input.t_max; input.t_min; input.endtime]*3600;
input.Qv = [input.Qv_max ; input.Qv_max ; input.Qv_min; input.Qv_min];
input.Qv = [input.tv,input.Qv];
end

% Numerics and Discretization
input.dz = input.H/(input.m-1);
input.z = (0:input.dz:input.H)';

% Noise
input.noisetime = (0:input.dt:input.endtime*3600)';
input.noiserand = input.noisesd*randn(length(input.noisetime),1)+input.noisemean;
input.noise = [input.noiserand, input.noisetime];

% Export
outputs = input;

return