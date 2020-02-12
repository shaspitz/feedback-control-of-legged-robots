clc
close all
clear all
%% Simulate Two Link Walker - Problem 1(d)
x0 = [0.2065;
    0.4130;
    -0.2052;
    -0.0172];

% init data vectors
t_vec = [];
x_vec = [];

% Loop for 10 steps
for i = 1:10
% Define time range to simulate the system
Tspan = [0 15] ;
t0 = 0 ; % Initial Time

% Initialize vectors
t_ode = []; x_ode = [];

% Define the event functions (stop integration when impact happens)
options = odeset('Events', @two_link_event) ;

% Simulate the system for each step
[t_ode, x_ode] = ode45(@two_link_dynamics, t0+Tspan, x0, options) ;

% Save simulation data
t_vec = [t_vec; t_ode] ;
x_vec = [x_vec; x_ode] ;

% Store row numbers of new steps for animation
t_I(i) = size(t_vec,1);

% Initialize xo and t for next step
x0(1:4) = two_link_impactdynamics(x_ode(end,1:4)); t0 = t_vec(end);
end

% Animate
animate_two_link_walker(t_vec, x_vec, 0.01, t_I)
