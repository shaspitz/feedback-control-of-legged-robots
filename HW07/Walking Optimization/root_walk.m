%% Feedback Control of Legged Robots - Homework #7
% Shawn Marshall-Spitzbart
% UC Berkeley

clear all; close all; clc; % starting with a clean slate

addpath('C:\Users\shawn\OneDrive\Documents\Berkeley\ME292B\HW07\gen')

%% Problem 1

% Initial seed for optimization
alpha0 = pi/6; alpha1 = 0; alpha2 = 0; alpha3 = 0;
beta0 = 0; beta1 = 0; beta2 = 0; beta3 = 0; 
x0 = [-0.3827;
    0.9239;
    2.2253;
    3.0107;
    0.5236;
    0.8653;
    0.3584;
    -1.0957;
    -2.3078;
    2.0323];
x_optim0 = [x0(1:10); alpha0; alpha1; alpha2; alpha3; ...
    beta0; beta1; beta2; beta3];

% Only linear constraint is to start 'x' position of hip at 0
Aineq = []; Bineq = [];
Aeq = [1 zeros(1,17)]; Beq = [0]; 
LB = [];
UB = [];

options = optimset('MaxFunEvals',20000,'MaxIter',20000);

% Set Desired Velocity
param.vd = 0.5;

[x_optim,obj_optim,feas_flag] = fmincon(@walk_obj,x_optim0,Aineq,Bineq,Aeq,Beq,LB,UB,@walk_cons,options,param);

% Display constraints
[cineq,ceq] = walk_cons(x_optim,param)

%% Simulate and Animate system

% Simulate Function for 5 steps now (alternative simulation function used)
[t_sim, x_sim] = sim_walk_alt(x_optim, param);

% Animate
animateThreeLink(t_sim, x_sim)

% Plots
plotting(t_sim,x_sim,x_optim);

%% Problem 2

% Report the fastest feasible periodic walking gait, and the
% slowest feasible periodic walking gait. Then plot the same
% information as above. 

% Set Desired Starting Velocity
param.vd = 0.7;

% Find fastest feasible gait by iterating through increasing desired velocities
% until optimization problem is infeasible
while feas_flag > 0
x_optim_fastest = x_optim; 
[x_optim,obj_optim,feas_flag] = fmincon(@walk_obj,x_optim0,Aineq,Bineq,Aeq,Beq,LB,UB,@walk_cons,options,param);
param.vd = param.vd + 0.1;
end

% Print fastest gait, its optimum values and constraints
fastest_gait_velocity = param.vd
x_optim_fastest
[cineq,ceq] = walk_cons(x_optim_fastest,param)

% Simulate Function for 5 steps now (alternative simulation function used)
[t_sim, x_sim] = sim_walk_alt(x_optim, param);

% Plots
plotting(t_sim,x_sim,x_optim);

%%
% Reset Desired Starting Velocity and feas flag
param.vd = 0.7;
feas_flag = 1;

% Find slowest gait by iterating through decreasing desired velocities
% until optimization problem is infeasible
while feas_flag > 0
x_optim_slowest = x_optim; 
[x_optim,obj_optim,feas_flag] = fmincon(@walk_obj,x_optim0,Aineq,Bineq,Aeq,Beq,LB,UB,@walk_cons,options,param);
param.vd = param.vd - 0.05;
end

% Print slowest gait, its optimum values and constraints
slowest_gait_velocity = param.vd
x_optim_slowest
[cineq,ceq] = walk_cons(x_optim_slowest,param)

% Simulate Function for 5 steps now (alternative simulation function used)
[t_sim, x_sim] = sim_walk_alt(x_optim, param);

% Plots
plotting(t_sim,x_sim,x_optim);

%%
% Placeholder for MATLAB plotting
placeholder = 1

%% Functions

function obj_sum = walk_obj(x,param)

    % Simulate walk
    [t_vec, x_vec] = sim_walk(x, param) ;
    
    alpha = x(11:14)';
    beta = x(15:18)';
    
    u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
    beta) + v_gen(s, alpha, beta));
    
    % cumulative cost
    x_TI = x_vec(end,1);
    
    % Implement the integral cost function 
    obj_sum = 0;
    for i = 1:length(t_vec)-1
    obj_sum = obj_sum + (1 / x_TI ) * norm( u(x_vec(i,:)', alpha,...
        beta) ) * (t_vec(i+1) - t_vec(i));
    end
    
end

function [cineq,ceq] = walk_cons(x,param)

    % unraveling the input variable
    x0 = x(1:10);
    alpha = x(11:14)';
    beta = x(15:18)';
    
    u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
        beta) + v_gen(s, alpha, beta));

    % Simulate 1 step to produce (non-decision variable) paramters that will 
    % feed into nonlinear constraints
    [t_vec, x_vec] = sim_walk(x, param);
    
    u_vec = zeros(size(x_vec, 1), 2);
    for row = 1:size(x_vec, 1)
        u_vec(row,1:2) = u(x_vec(row,:)', alpha, beta);
    end
    
    for row = 1:size(x_vec, 1)
        Fst(row,:) = Fst_gen(x_vec(row,:)', u_vec(row,:)');
    end

    % inequality constraints (nonlinear)
    cineq = [];
    
    % constraint (a) - Unilateral ground constraints
    cineq = [cineq; max(-Fst(:, 2))];
    
    % constraint (b) - Friction cone constraints
    mu_s = 0.75;
    cineq = [cineq;  max( abs(Fst(:, 1)./Fst(:, 2))) - mu_s ];
    
    % constraint (c) - Desired Average Speed Constraints
    x_TI = x_vec(end,1);
    TI = t_vec(end); 
    cineq = [cineq; -x_TI / TI + param.vd];
    
    % equality constraints
    ceq = [];
    
    % constraint (d) - Periodicity constraints 
    R = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 1 0;
    0 0 1 0 0;
    0 0 0 0 1];
    ceq = [ceq; x0(2:5,1) - R(2:end,2:end) * x_vec(end,2:5)']; % Dont include 'x' since robot moves forward after step
    ceq = [ceq; x0(6:10,1) - R * dqPlus_gen(x_vec(end,:)')];
    
end

function [t_ode, x_ode] = sim_walk(x, param)

% unraveling the input variable
q0 = x(1:10);
alpha = x(11:14)';
beta = x(15:18)';

% Define u and ds
u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
    beta) + v_gen(s, alpha, beta));

% Define state function to integrate
ds = @(t,s,alpha,beta) f_gen(s) + g_gen(s) * u(s, alpha, beta);

% Define time range to simulate the system
tspan = [0 10] ;

% Define the event functions (stop integration when impact happens)
options = odeset('Events', @three_link_event);

% Simulate the system for each step
[t_ode,x_ode] = ode45(@(t,s) ds(t,s,alpha,beta),tspan,q0,options);

end
function [value,isterminal,direction] = three_link_event(t,x)

value = x(3) + x(5) - pi - pi/8; % detect when phi - 2*theta == 0 (approx)
isterminal = 1 ; % stop integration when value == 0
direction = 1 ; % detect zero when function is increasing

end

function [t_vec, x_vec] = sim_walk_alt(x, param)

% unraveling the input variable
q0 = x(1:10);
alpha = x(11:14)';
beta = x(15:18)';

% Define u and ds
u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
    beta) + v_gen(s, alpha, beta));

% Define state function to integrate
ds = @(t,s,alpha,beta) f_gen(s) + g_gen(s) * u(s, alpha, beta);

% Initialize vectors
t_vec = []; x_vec = [];

t0 = 0 ; % Initial Time

% Impact map
R = [1 0 0 0 0;
0 1 0 0 0;
0 0 0 1 0;
0 0 1 0 0;
0 0 0 0 1];

    % Loop for 5 steps
    for i = 1:5

    % Define time range to simulate the system
    tspan = [0 10] ;

    % Define the event functions (stop integration when impact happens)
    options = odeset('Events', @three_link_event);

    % Simulate the system for each step
    [t_ode,x_ode] = ode45(@(t,s) ds(t,s,alpha,beta),t0+tspan,q0,options);

    % Save simulation data
    t_vec = [t_vec; t_ode] ;
    x_vec = [x_vec; x_ode] ;

    % Initialize xo and t for next step
    q0(1:5,1) = R * x_ode(end,1:5)';
    q0(6:10,1) = R * dqPlus_gen(x_ode(end,:)');
    t0 = t_vec(end);

    end

end

function plotting(t_sim,x_sim,x_optim)

% Theta1 vs dTheta1

% Transformation matrix:
T = [1 0 0 0 0;
     0 1 0 0 0;
     0 0 1 0 1;
     0 0 0 1 1;
     0 0 0 0 1];
d = [0;
     0;
     -pi;
     -pi;
     0];

% First convert data to theta coordinates
q_tild = zeros(size(x_sim, 1), 5);
dq_tild = zeros(size(x_sim, 1), 5);

for row = 1:size(x_sim, 1)
    q_tild(row, 1:5) = T * x_sim(row,1:5)' + d;
    dq_tild(row, 1:5) = T * x_sim(row,6:10)' + zeros(5,1);
end

% Plot
figure()
plot(dq_tild(:, 3) , q_tild(:, 3))
title('Theta1 vs dTheta1')
xlabel('dTheta1')
ylabel('Theta1')

% u1 and u2 vs time
u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
beta) + v_gen(s, alpha, beta));

uplot = zeros(size(x_sim, 1), 2);
for row = 1:size(x_sim, 1)
    uplot(row,1:2) = u(x_sim(row,1:10)', x_optim(11:14)', x_optim(15:18)');
end

u1plot = uplot(:,1);
u2plot = uplot(:,2);

% Plot
figure()
plot(t_sim(:,1), u1plot(:,1))
hold on
plot(t_sim(:,1), u2plot(:,1))
legend('u1', 'u2', 'Location', 'Best')
title('u1 and u2 vs Time')
xlabel('Time')
ylabel('u1 and u2')

% Fst vs time
Fst_plot = zeros(size(x_sim, 1), 2);
for row = 1:size(x_sim, 1)
    Fst_plot(row,1:2) = Fst_gen(x_sim(row,:)', uplot(row,:)');
end

Fst_plot_tot = sqrt(Fst_plot(:,1).^2 + Fst_plot(:,2).^2);

% Plot
figure()
plot(t_sim(:,1), Fst_plot_tot)
hold on
plot(t_sim(:,1), Fst_plot(:,1))
plot(t_sim(:,1), Fst_plot(:,2))
legend('Fst Magnitude', 'Fst (horizontal component)', 'Fst (vertical component)','Location', 'Best')
title('Fst vs Time')
xlabel('Time')
ylabel('Fst')

end

