clc;
clear
%% ============================= Parameters ===============================
% Snake parameters (mass: m; length: 2l)
N = [7;15;50];

for ii=1:3
    clearvars -except N ii 
    param.N  = N(ii);
    param.m  = 0.406;
    param.l  = 0.0525;
    param.g  = 9.81;
    param.priemer = 0.30;
    param.priemerInfluence = 0.10;
    dlzkaPotrubia = 3;
    param.d = param.priemer - 2*param.l;
    pocetSnimkov = 50;
    param.dt = 0.01;
    % Friction coefficients:
    param.ct = 0.015;%0.02
    param.cn = 0.03;
    param.ut = 0.15; %zatial 0.15 ku 0.3 najlepsie
    param.un = 0.25;
    param.ctPipe = 0.015;
    param.utPipe = 0.2;
    param.umax = 3;
    param.qmax = 400*param.dt; % 400 st/s ale v case dt to je 400*dt 
    param.Erub = 400000; %last 400000 %95000
    param.vrub = 0.49; 
    param.tlmic = .5;
    param.pruzina = 5;
    param.trenie = 1;   % 0 - Coulomb, 1 - viscous
    param.kontakt = 1;   % 0 - bez, 1 - s
    param.minLinkVel = 0.001;   % [mm]
    param.dimensionPlot3D = 0;  % 0 - 2D, 1 - 3D
    param.resultsShow = 0;  % 0 - show simulation, 1 - show graphs 
    
    
    % Optimization
    optimization = 0;   % 0 - without, 1 - with
    
    % Controler parameters:
    param.kp  = 25; %%%***
    param.kd  = 10; %%%***
    
    % Simulation time:
    t=0:param.dt:60;
    
    % Reference trajectory parameters:
    
    % 0.3981    0.6936    0.4914
    param.alfa = 0.3981;%deg2rad(56.4012);     % 52.11
    param.omega = 0.6936;%deg2rad(49.7791);    % 0.13
    param.delta = 0.4914;%deg2rad(1.4074);  % 0.5041
    param.offset = 0;
    
    % Trajectory
    for i=1:param.N-1
        fi = param.alfa*sin((param.omega*t+(i-1)*param.delta));
        fi_required{i} = fi;
        fi_reference(i) = fi(1);
    end
    
    % Initial values
    theta       = zeros(param.N,1);
    thetaDot    = zeros(param.N,1);
    fi          = zeros(1,param.N-1); % teraz ide rovno od zaciatku, inak toto: fi_reference
    fiDot       = zeros(param.N-1,1);
    p           = zeros(2,1);
    pDot        = zeros(2,1);
    
    qa          = fi';
    qu          = [theta(param.N);p(1);p(2)];
    qaDot       = fiDot;
    quDot       = [thetaDot(param.N);pDot(1);pDot(2)];
    x0          = [qa;qu;qaDot;quDot;zeros(2*param.N,1)];%;zeros(2*param.N,1)
    %x0          = [qa;qu;qaDot;quDot];
    
    %% Solve
    if optimization == 0
        [T,X] = ode45(@(t,y)dynamicModel_last(t,y,param),t,x0);
        for k=1:2*param.N
           % contactForces(:,k) = diff(X(:,2*param.N+k));
        end
        %{
        for k=1:param.N
            %torque(:,k) = diff(X(:,param.N+k));
            torque(:,k) = X(:,param.N+k);
        end
        %}
        subplot(3,1,ii)
        run('animacia.m')
        hold on
    end
    
end

