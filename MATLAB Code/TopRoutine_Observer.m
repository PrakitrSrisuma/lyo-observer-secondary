% ==============================================================================
% This is a top-level routine for state estimation.
%
% Created by Prakitr Srisuma, 
% PhD, Braatz Group (ChemE) & 3D Optical Systems Group (MechE), MIT.
% ==============================================================================
close all; clear; clc;

%% Pre-simulation
% Add paths
addpath('Input Data', 'Validation Data', 'Saved Data', 'PDEs', 'Events', 'Calculations');

% Define and extract important input data
input = get_input_data; 
ip = input_processing(input);
m = ip.m;

% Mode - set on or off
Obs = 'on';  % simple example
ObsDes0 = 'off';  % sim-based design; Fig. 4
ObsDes1 = 'off';  % sim-based design; Figs. 5, 6
ObsDes2 = 'off';  % sim-based design; Fig. 7
ObsDes3 = 'off';  % sim-based design with noise; Fig. 8
ObsExp = 'off';  % exp-based design; Fig. 9
ObsExp2 = 'off';  % exp-based design; Fig. 10  
ObsMW = 'off';  % for microwave lyo; Fig. 11


%% State Estimation
switch Obs
case 'on'

% Observer gain
input = get_input_data; 
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;1.1*T0;zeros(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
ec = abs(cs_est - cs);
eT = abs(T_est - T);
tcon = t(find(ec<0.02*ec(1),1));

% Plot
figure; plot(t,T(:,end),'-.r','linewidth',4); hold on; plot(t,T_est(:,end),'-b','linewidth',2)
legend('True state','Estimated state'); ylabel('Bound water concentration (kg water/kg solid)'); xlabel('Time (h)')
figure; plot(t,cs(:,end),'-.r','linewidth',4); hold on; plot(t,cs_est(:,end),'-b','linewidth',2)
legend('True state','Estimated state'); ylabel('Bound water concentration (kg water/kg solid)'); xlabel('Time (h)')

end


%% Simulation-based Observer Design
switch ObsDes0
case 'on'

Obs0 = figure;

% Observer gain analysis
input = get_input_data; 
npoint = 20;
LT = -logspace(-6,-2,npoint);
Lc = logspace(-8,-6,npoint);
[X,Y] = meshgrid(LT,Lc);
nrun = length(LT)*length(Lc);
t0 = [];

for i = 1:length(LT)
    for j = 1:length(Lc)
        input.LT = LT(i);
        input.Lc = Lc(j);
        ip = input_processing(input);
    
        % Initial conditions
        T0 = ip.T0*ones(m,1);
        cs0 = ip.cs0*ones(m,1);
        
        % ODE solver setup
        tspan = (0:input.dt:input.endtime*3600)';
        y0 = [T0;cs0;T0;0.0314*ones(m,1)];
        opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);
        
        % Simulation
        tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;
        
        % Post-processing
        t = t/3600;
        T = y(:,1:m);
        cs = y(:,m+1:2*m);
        T_est = y(:,2*m+1:3*m);
        cs_est = y(:,3*m+1:end);
        Tavg = mean(T,2);
        cs_avg = mean(cs,2);
        Tavg_est = mean(T_est,2);
        cs_avg_est = mean(cs_est,2);
        error = abs(cs_avg-cs_avg_est);
        t0(i,j) = t(find(error<.02*error(1),1));

    end
end

surf(X,Y,t0','FaceAlpha',0.9,'edgecolor','interp','FaceColor','flat')
h=gca;
graphics_setup('1by1.5')
set(h,'xscale','log','yscale','log')
xlabel('{\it L_T}'); ylabel('{\it L_c}'); zlabel('Convergence time (h)')
view(125,25); colormap(winter)
cb = colorbar('location','eastoutside'); 
set(cb,'YTick',(0:2:10),'ylim',[0 10])

end


%% Simulation-based Observer Design
switch ObsDes1
case 'on'

% For average values
input = get_input_data;
input.LT = -1e-6;
input.Lc = 5e-7;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;1.1*T0;0.0314*ones(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);


% For spatiotemporal properties
input = get_input_data;  % default inputs
input.Lc = 5e-7;
input.LT = 1e-6;
ip = input_processing(input);

% ODE solver setup
tspan = (0:1500:input.endtime*3600)';
y0 = [T0;cs0;1.1*T0;0.6415*ones(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t2,y2] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t2 = t2/3600;
T2 = y2(:,1:m);
cs2 = y2(:,m+1:2*m);
T_est2 = y2(:,2*m+1:3*m);
cs_est2 = y2(:,3*m+1:end);
Tavg2 = mean(T2,2);
cs_avg2 = mean(cs2,2);
Tavg_est2 = mean(T_est2,2);
cs_avg_est2 = mean(cs_est2,2);
eT = abs(T_est2-T2);
ec = abs(cs_est2-cs2);

% Plot
figure
Obs1T = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
nexttile
plot(t,Tavg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,Tavg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','best'); ylabel('Average temperature (K)'); xlabel('Time (h)')
h.ItemTokenSize(1) = 15;
xlim([0 12])
xticks(0:2:12)
graphics_setup('1by2b')
set(gca,'XMinorTick','on','YMinorTick','on')
text(0.02,.95,'(A)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');

nexttile;
surf(ip.z*100, t2, eT,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error (K)'})
ylim([0 12])
yticks(0:4:12)
graphics_setup('1by2b')
set(gca, 'ydir', 'reverse','XMinorTick','on','YMinorTick','on')
colormap(copper)
text(0,.95,'(B)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');

figure
Obs1c = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
nexttile
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_est,'-b','linewidth',2); 
h = legend('True state','Estimated state','location','best'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 15;
xlim([0 12])
xticks(0:2:12)
graphics_setup('1by2')
text(0.02,.95,'(A)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');

nexttile;
surf(ip.z*100, t2, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
graphics_setup('1by2')
ylim([0 12])
yticks(0:4:12)
set(gca, 'ydir', 'reverse')
colormap(copper)
text(0,.95,'(B)','Units','normalized','FontSize', 9 ,'fontweight', 'bold');

end


%% Simulation-based Observer Design - Various Conditions
switch ObsDes2
case 'on'

figure; Obs2 = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

%%% Case A_min
ObsA = tiledlayout(Obs2,1,2);
ObsA.Layout.Tile = 1; 
input = get_input_data; 
input.A = 7.8e-5;
input.Ea = 0;
input.endtime = 15;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsA)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
xticks(0:2:10)
xlim([0 10])
set(gca,'XMinorTick','on','YMinorTick','on')  
graphics_setup('3by4')
title({'';'';'';''})
text(-0.35,1.3,'(A) Variation in frequency factor','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
text(-0.35,1.17,'1. Low case, {\it A} = 7.8\times10^{–5} s^{-1}','Units','normalized','FontSize', 8 ,'fontweight', 'bold');

nexttile(ObsA)
tspan = (0:1500:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 10])
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)


%%% Case A_max
ObsA = tiledlayout(Obs2,1,2);
ObsA.Layout.Tile = 2; 
input = get_input_data;
input.A = 1.1e-4;
input.Ea = 0;
input.endtime = 15;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsA)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
xticks(0:2:8)
xlim([0 8])
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('3by4')
text(-0.35,1.17,'2. High case, {\it A} = 1.1\times10^{–4} s^{-1}','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
nexttile(ObsA)
tspan = (0:1500:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 8])
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)


%%% Case B_min
ObsB = tiledlayout(Obs2,1,2);
ObsB.Layout.Tile = 3; 
input = get_input_data;
input.Ea = 5920;
input.endtime = 15;
ip = input_processing(input);


% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsB)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
xticks(0:1:4)
xlim([0 4])
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('3by4')
title({'';'';'';'';''})
text(-0.35,1.3,'(B) Variation in activation energy','Units','normalized','FontSize', 8,'fontweight', 'bold');
text(-0.35,1.15,'1. Low case, {\it E_a} = 5,920 J/mol','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
nexttile(ObsB)
tspan = (0:1000:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 4])
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)


%%% Case B_max
ObsB = tiledlayout(Obs2,1,2);
ObsB.Layout.Tile = 4; 
input = get_input_data;
input.Ea = 13416;
input.endtime = 100;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsB)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('3by4')
xlim([0 50])
text(-0.35,1.15,'2. High case, {\it E_a} = 13,416 J/mol','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
nexttile(ObsB)
tspan = (0:8000:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 50])
yticks(0:25:50)
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)


%%% Case C_min
ObsC = tiledlayout(Obs2,1,2);
ObsC.Layout.Tile = 5; 
input = get_input_data; 
input.cs0 = 0.0314;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsC)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
xticks(0:1:4)
xlim([0 4])
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('3by4')
title({'';'';'';''})
text(-0.35,1.3,'(C) Variation in initial concentration','Units','normalized','FontSize', 8,'fontweight', 'bold');
text(-0.35,1.15,'1. Low case, {\itc_s}_{,0} = 0.0314 kg water/kg solid','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
nexttile(ObsC)
tspan = (0:1000:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 4])
zlim([0 0.2])
clim([0 .2])
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)


%%% Case C_max
ObsC = tiledlayout(Obs2,1,2);
ObsC.Layout.Tile = 6; 
input = get_input_data;
input.cs0 = 0.6415;
input.endtime = 100;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);

% Plot
nexttile(ObsC)
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state','location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
graphics_setup('3by4')
xticks(0:2:10)
xlim([0 10])
set(gca,'XMinorTick','on','YMinorTick','on') 
text(-0.35,1.15,'2. High case, {\itc_s}_{,0} = 0.6415 kg water/kg solid','Units','normalized','FontSize', 8 ,'fontweight', 'bold');
nexttile(ObsC)
tspan = (0:2000:input.endtime*3600)';
tic; [t,y] = ode15s (@(t,y) PDE_Observer(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 10])
set(gca, 'ydir', 'reverse')
colormap(copper)
graphics_setup('3by4')
view(-45,20)

end


%% Simulation-based Observer Design - Noise
switch ObsDes3
case 'on'

figure; Obs3 = tiledlayout(2,3,'TileSpacing','loose','Padding','compact');

% Input
input = get_input_data; 
input.Lc = 5e-7;
ip = input_processing(input);
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);
T_noisy = load('Sim_T_noisy').T_noisy;
Tavg_noisy = mean(T_noisy,2);
noise = load('Sim_noise').noise;

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer_Noise(t,y,T_noisy,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
T_measured = mean(T_noisy(:,1:end-1),2);

% Plot
nexttile
title({'';''})
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(.45,.65,'{\it L_c}  = 5\times10^{–7}','Units','normalized','FontSize', 7);
text(-0.35,1.10,'(A) Observer gain selection for noisy temperature measurement','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
xlim([0 12])
xticks(0:2:12)
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('2by3')


%%% Simulation 2
input = get_input_data;  % default inputs
input.Lc = 2e-7;
ip = input_processing(input);

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer_Noise(t,y,T_noisy,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);

% Plot
nexttile
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
h = legend('True state','Estimated state'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(.45,.65,'{\it L_c}  = 2\times10^{–7}','Units','normalized','FontSize', 7);
xlim([0 12])
xticks(0:2:12)
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('2by3')


%%% Simulation 3
input = get_input_data; 
input.Lc = 1e-7;
ip = input_processing(input);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer_Noise(t,y,T_noisy,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
T_measured = mean(T_noisy(:,1:end-1),2);

% Plot
nexttile
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
title({'';''})
h = legend('True state','Estimated state'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(.45,.65,'{\it L_c}  = 1\times10^{–7}','Units','normalized','FontSize', 7);
xlim([0 12])
xticks(0:2:12)
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('2by3')


%%% Simulation 3
input = get_input_data; 
ip = input_processing(input);
tau_0 = cal_timeconstant(ip).tau_0;

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer_Scheduling(t,y,T_noisy,tau_0,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
T_measured = mean(T_noisy(:,1:end-1),2);

% Plot
nexttile
plot(t,cs_avg,'linewidth',8,'Color',[1, 0, .7, 0.4]); hold on; plot(t,cs_avg_est,'-b','linewidth',2)
title({'';''})
h = legend('True state','Estimated state'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(-0.35,1.11,'(B) Gain scheduling: {\itL_c} changed from 5\times10^{–7} to 1\times10^{–7} at {\itt} = 4{\it𝝉}','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
xlim([0 12])
xticks(0:2:12)
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('2by3')


input = get_input_data; 
ip = input_processing(input);
tspan = (0:1500:input.endtime*3600)';

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_Observer_Scheduling(t,y,T_noisy,tau_0,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
Tavg = mean(T,2);
cs_avg = mean(cs,2);
Tavg_est = mean(T_est,2);
cs_avg_est = mean(cs_est,2);
T_measured = mean(T_noisy(:,1:end-1),2);

% Plot
nexttile
eT = abs(T_est-T);
ec = abs(cs_est-cs);
surf(ip.z*100, t, ec,'linewidth',1,'FaceColor','white','EdgeColor','interp','MeshStyle','row'); xlabel('Position (cm)'); ylabel('Time (h)'); zlabel({'Estimation error';'(kg water/kg solid)'})
ylim([0 12])
yticks(0:4:12)
set(gca, 'ydir', 'reverse')
graphics_setup('2by3')

colormap(copper)

nexttile
plot(T_noisy(:,end)/3600,T_measured,'-k'); 
xlim([0 12])
xticks(0:2:12)
ylabel('Average temperature (K)'); xlabel('Time (h)')
set(gca,'XMinorTick','on','YMinorTick','on') 
text(-0.24,1.12,'(C) Noisy temperature measurement','Units','normalized','FontSize', 8 ,'fontweight', 'bold' );
graphics_setup('2by3')

end


%% Experimental-based Observer Design
switch ObsExp
case 'on'

% Validation data 1
t_sim = load('Data1Sim_t.mat').t;
T_sim = load('Data1Sim_Temp.mat').T;
data = load('Data1.mat').cs;
cs_exp = data(:,2);
t_exp = data(:,1);
input = get_input_data_exp1; 
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (ip.cs0-0.5*ip.cs0)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i,:)',ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
ObsExp = tiledlayout(2,2,'TileSpacing','loose','Padding','compact');
nexttile
plot(t_sim,cs_est,'-b','linewidth',2);
hold on; plot(t_exp,cs_exp,'sm','MarkerSize',4,'MarkerFaceColor','m');
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(0.05,0.1,'(A)','Units','normalized','FontSize', 11 ,'fontweight', 'bold' );
ylim([0 .8])
graphics_setup('2by2')
set(gca,'XMinorTick','on','YMinorTick','on') 

%============================

% Validation data 5
t_sim = load('Data5Sim_t.mat').t;
T_sim = load('Data5Sim_Temp.mat').T;
data = load('Data5.mat').cs;
cs_exp = data(:,2);
t_exp = data(:,1);
input = get_input_data_exp5; 
% input.LT = -3e-5;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (.27)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i,:)',ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
nexttile
plot(t_sim,cs_est,'-b','linewidth',2);
hold on; plot(t_exp,cs_exp,'sm','MarkerSize',4,'MarkerFaceColor','m');
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(0.05,0.1,'(B)','Units','normalized','FontSize', 11 ,'fontweight', 'bold' );
graphics_setup('2by2')
xlim([0 22])
set(gca,'XMinorTick','on','YMinorTick','on') 

%============================

% Validation data 4
cs = []; cs_est = []; T = [];  T_est = [];
t_sim = load('Data4Sim_t.mat').t;
T_sim = load('Data4Sim_TB.mat').TB;
cs_exp = load('Data4_cs.mat').cs;
cs_avg = load('Data4Sim_csavg.mat').cs_avg;
input = get_input_data_exp4; 
input.obs = 'OneTemp';
input.LT = -5e-3;
input.Lc = 1e-4;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (.06)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
ts = 0;
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i),ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    ts = [ts;t(end)];
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
ts = ts/3600;
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
nexttile
plot(t_sim,cs_est,'-b','linewidth',2);
hold on
errorbar(cs_exp(:,1),cs_exp(:,2),cs_exp(:,3),'sm','MarkerSize',4,'MarkerFaceColor','m','CapSize',6);
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
ylim([0 inf])
graphics_setup('2by2')
xticks(0:1:7)
xlim([0 7])
text(0.05,0.1,'(C)','Units','normalized','FontSize', 11 ,'fontweight', 'bold' );
set(gca,'XMinorTick','on','YMinorTick','on') 

%============================

% Validation data 3
cs = []; cs_est = []; T = [];  T_est = [];
t_sim = load('Data3Sim_t.mat').t;
T_sim = load('Data3Sim_TB.mat').TB;
cs_exp = load('Data3_cs.mat').cs;
cs_avg = load('Data3Sim_csavg.mat').cs_avg;
input = get_input_data_exp3; 
input.obs = 'OneTemp';
input.LT = -5e-3;
input.Lc = 1e-4;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (.0314)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
tend = 9;  % hours
ts = 0;
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i),ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    ts = [ts;t(end)];
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
ts = ts/3600;
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
nexttile
plot(t_sim,cs_est,'-b','linewidth',2); 
hold on;  
errorbar(cs_exp(:,1),cs_exp(:,2),cs_exp(:,3),'sm','MarkerSize',4,'MarkerFaceColor','m','CapSize',6);
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
graphics_setup('2by2')
xticks(0:1:6)
ylim([0 inf])
text(0.05,0.1,'(D)','Units','normalized','FontSize', 11 ,'fontweight', 'bold' );
set(gca,'XMinorTick','on','YMinorTick','on') 

end


%% Experimental-based Observer Design
switch ObsExp2
case 'on'

% Validation data 1
t_sim = load('Data1Sim_t.mat').t;
T_sim = load('Data1Sim_Temp.mat').T;
data = load('Data1.mat').cs;
cs_exp = data(:,2);
t_exp = data(:,1);
input = get_input_data_exp1; 
input.h = input.h-0.1*input.h;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (ip.cs0-0.5*ip.cs0)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i,:)',ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
ObsExp2 = tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
nexttile
plot(t_sim,cs_est,'-b','linewidth',2);
hold on; plot(t_exp,cs_exp,'sm','MarkerSize',4,'MarkerFaceColor','m');
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
ylim([0 .8])
graphics_setup('1by2')
set(gca,'XMinorTick','on','YMinorTick','on') 

%============================

% Validation data 5
t_sim = load('Data5Sim_t.mat').t;
T_sim = load('Data5Sim_Temp.mat').T;
data = load('Data5.mat').cs;
cs_exp = data(:,2);
t_exp = data(:,1);
input = get_input_data_exp5; 
input.h = input.h + 0.1*input.h;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = (.27)*ones(m,1);
y0 = [T0;cs0];

% ODE solver setup
tspan = t_sim*3600;
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
for i = 1:length(t_sim)-1
    [t,y] = ode15s (@(t,y) PDE_Observer_Exp(t,y,T_sim(i,:)',ip), [tspan(i) tspan(i+1)], y0, opts_ode);
    T(i,:) = y(end,1:m);
    cs(i,:) = y(end,m+1:2*m);
    y0 = [T(i,:),cs(i,:)];   
end

% Post-processing
cs = [cs0';cs];
cs_est = mean(cs,2);
T = [T0';T];
T_est = mean(T,2);

% Plot
nexttile
plot(t_sim,cs_est,'-b','linewidth',2);
hold on; plot(t_exp,cs_exp,'sm','MarkerSize',4,'MarkerFaceColor','m');
h = legend('Estimated state', 'True state', 'location','northeast'); ylabel({'Average concentration';'(kg water/kg solid)'}); xlabel('Time (h)')
h.ItemTokenSize(1) = 10;
text(0.05,0.1,'(B)','Units','normalized','FontSize', 11 ,'fontweight', 'bold' );
graphics_setup('2by2')
xlim([0 22])
set(gca,'XMinorTick','on','YMinorTick','on') 
graphics_setup('1by2')

end


%% State Estimation
switch ObsMW
case 'on'

% Observer gain
input = get_input_data;
input.endtime = 8;
input.mode = 'MFD';
input.K = 1000;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);
noise = ip.noise;

% ODE solver setup
tspan = (0:input.dt:input.endtime*3600)';
y0 = [T0;cs0;T0;0.0314*ones(m,1)];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t,y] = ode15s (@(t,y) PDE_ObserverFB(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t = t/3600;
T = y(:,1:m);
cs = y(:,m+1:2*m);
Tavg = mean(T,2);
T_est = y(:,2*m+1:3*m);
cs_est = y(:,3*m+1:end);
cs_avg = mean(cs,2);


%%% CFD
input = get_input_data; 
input.mode = 'CFD';
input.endtime = 8;
ip = input_processing(input);

% Initial conditions
T0 = ip.T0*ones(m,1);
cs0 = ip.cs0*ones(m,1);

% ODE solver setup
tspan = (0:ip.dt:ip.endtime*3600)';
y0 = [T0;cs0];
opts_ode = odeset('RelTol',1e-10,'AbsTol',1e-10);
opts_terminate = odeset('Event',@(t,y) event_desorption_completes(t,y,input),'RelTol',1e-10,'AbsTol',1e-10);

% Simulation
tic; [t2,y2] = ode15s (@(t,y) PDE_ModelFVM(t,y,ip), tspan, y0, opts_ode); toc;

% Post-processing
t2 = t2/3600;
T2 = y2(:,1:m);
Tavg2 = mean(T2,2);
cs2 = y2(:,m+1:2*m);
cs_avg2 = mean(cs2,2);

% Plotting
ObsMFD = tiledlayout(1,2,'TileSpacing','loose','Padding','compact'); 
nexttile
plot(t2,cs_avg2,'linewidth',6,'Color',[.2, 0.7, .7, 0.4]); hold on; plot(t,cs_avg,'linewidth',6,'Color',[1, 0, .7, 0.4]); hold on; 
plot(t,cs_est(:,end),'-b','linewidth',2); 
h = legend('True state (conventional)','True state (microwave)','Estimated state','location','best'); ylabel({'Average concentration'; '(kg water/kg solid)'}); xlabel('Time (h)')
text(.8,.35,'(A)','Units','normalized','FontSize', 11,'fontweight', 'bold');
graphics_setup('1by2')
set(gca,'XMinorTick','on','YMinorTick','on') 
nexttile
plot(t,Tavg,'-r','linewidth',2)
ylabel('Average temperature (K)'); xlabel('Time (h)')
text(.8,.35,'(B)','Units','normalized','FontSize', 11,'fontweight', 'bold');
h.ItemTokenSize(1) = 10;
graphics_setup('1by2')
set(gca,'XMinorTick','on','YMinorTick','on') 

end
