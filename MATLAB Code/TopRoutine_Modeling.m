% ==============================================================================
% This is a top-level routine for model validation.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Validation Data', 'Saved Data', 'PDEs', 'Events', 'Calculations');

% Define and extract important input data
input = get_input_data;  % default inputs
ip = input_processing(input);
m = ip.m;

% This corresponds to Fig. 3 in the paper.
figure
Exp = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');

%% Validation data 1
filename = 'Data1.mat';
cs_exp = load(filename).cs;
input = get_input_data_exp1;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tend = 9;  % hours
tspan = (0:ip.dt:tend*3600)';
y0 = [T0;cs0];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Solve the ODEs
tic; [t,y] = ode15s (@(t,y) PDE_ModelFVM(t,y,ip), tspan, y0, opts_ode); toc;
tic; [t_e,y_e] = ode15s (@(t,y) PDE_ModelFVM(t,y,ip), cs_exp(:,1)*3600, y0, opts_ode); toc;

% Output
t = t/3600;
cs = y(:,m+1:2*m);
cs_avg = mean(cs,2);
cs_e = mean(y_e(:,m+1:2*m),2);
error = cs_e - cs_exp(:,2);

% Plotting
nexttile
plot(t,cs_avg,'-b','linewidth',2); hold on; plot(cs_exp(:,1),cs_exp(:,2),'sm','MarkerSize',4,'MarkerFaceColor','m')
h = legend('Our model','Experimental data');
h.ItemTokenSize(1) = 15;
ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
xticks(0:2:10)
set(gca,'XMinorTick','on','YMinorTick','on') 
title({'';''})
text(-0.23,1.13,'(A) Validation with experimental data','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
graphics_setup('2by2')


%% Validation data 2
filename = 'Data2.mat';
cs_exp = load(filename).cs;
input = get_input_data_exp2; 
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tend = 12;  % hours
tspan = (0:ip.dt:tend*3600)';
y0 = [T0;cs0];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Solve the ODEs
tic; [t,y] = ode15s (@(t,y) PDE_ModelFVM(t,y,ip), tspan, y0, opts_ode); toc;

% Output
t = t/3600;
cs = y(:,m+1:2*m);
cs_avg = mean(cs,2);

% Plotting
nexttile
plot(t,cs_avg,'-b','linewidth',2); hold on; plot(cs_exp(:,1),cs_exp(:,2),'-om','linewidth',1,'MarkerSize',3,'MarkerFaceColor','m')
h = legend('Our model','High-fidelity model');
h.ItemTokenSize(1) = 15;
set(gca,'XMinorTick','on','YMinorTick','on') 
ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
title({'';''})
text(-0.23,1.13,'(B) Validation with high-fidelity model','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
graphics_setup('2by2')


%% Validation data 3
filename = 'Data3_cs.mat';
filename2 = 'Data3_T.mat';
cs_exp = load(filename).cs;
T_exp = load(filename2).T;
input = get_input_data_exp3;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);
pw0 = ip.pw0*ones(m,1);
pin0 = ip.pin0*ones(m,1);

% ODE solver setup
tend = 6;  % hours
tspan = (0:ip.dt:tend*3600)';
y0 = [T0;cs0];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Solve the ODEs
tic; [t,y] = ode15s (@(t,y) PDE_ModelFVM(t,y,ip), tspan, y0, opts_ode); toc;

% Output
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
Tavg = T(:,end);
TB = Tavg;
cs_avg = mean(cs,2);

% Plotting
nexttile; plot(t,Tavg,'-b','linewidth',2); hold on; plot(T_exp(:,1),T_exp(:,2),'sm','MarkerSize',4,'MarkerFaceColor','m')
h = legend('Our model','Experimental data','location','best');
h.ItemTokenSize(1) = 15;
set(gca,'XMinorTick','on','YMinorTick','on') 
ylabel({'Bottom temperature (K)'}); xlabel('Time (h)')
xticks(0:1:6)
title({'';''})
text(-0.23,1.13,'(C) Validation with experimental data: temperature and concentration','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
graphics_setup('2by2')
nexttile; plot(t,cs_avg,'-b','linewidth',2); hold on; errorbar(cs_exp(:,1),cs_exp(:,2),cs_exp(:,3),'sm','MarkerSize',4,'MarkerFaceColor','m','CapSize',6)
h = legend('Our model','Experimental data','location','best');
ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
ylim([0 inf])
xticks(0:1:6)
h.ItemTokenSize(1) = 15;
graphics_setup('2by2')
set(gca,'XMinorTick','on','YMinorTick','on') 

