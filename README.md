# Introduction
This repository contains software for the paper "Real-time estimation of bound water concentration during lyophilization with temperature-based state observers" published in "International Journal of Pharmaceutics". To use the code, see "Manual" in "MATLAB Code".

# Cite
If any material in this repository is useful for your work, please cite the following article:
https://doi.org/10.1016/j.ijpharm.2024.124693 or https://arxiv.org/abs/2407.13844

Srisuma, P., Barbastathis, G., Braatz, R.D., 2024. Real-time estimation of bound water concentration during lyophilization with temperature-based state observers. International Journal of Pharmaceutics, 124693. doi:10.1016/j.ijpharm.2024.124693.

# Example of Simulating a State Observer
```
input = get_input_data; 
ip = input_processing(input);
m = ip.m;
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;1.1*T0;zeros(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); 
```
