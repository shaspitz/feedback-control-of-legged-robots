clc
close all
clear all

%% ME 292B - Homework #3 Problem 1

%% Problem 1(d)

% Constants and input defined in stance_dynamics function

% Define time range to simulate the system
Tspan = linspace(0,10,100) ;
t0 = 0 ; % Initial Time

% Initial state vector of the system (ball position ; ball velocity)
x0 = [5;
      0] ;

% Array to store data for all bounces
t_vec = [] ; x_vec = [] ;

% Simulate 10 bounces
for j=1:10

    % Define the events functions (stop integration when stance stage reached)
    options1 = odeset('Events', @ground_contact) ;

    % Simulate the system in flight
    [t_ode x_ode] = ode45(@flight_dynamics, t0+Tspan, x0, options1) ;

    % Save simulation data
    t_vec = [t_vec; t_ode] ;
    x_vec = [x_vec; x_ode] ;

    % Initialize xo,t for stance stage
    x0 = [x_ode(length(x_ode),1); x_ode(length(x_ode),2)]; t0 = t_vec(end);

    % Define the events functions (stop integration when flight stage reached)
    options2 = odeset('Events', @leaves_ground) ;

    % Simulate the system in stance stage
    [t_ode x_ode] = ode45(@stance_dynamics, t0+Tspan, x0, options2) ;

    % Save simulation data
    t_vec = [t_vec; t_ode] ;
    x_vec = [x_vec; x_ode] ;

    % Initialize xo,t for flight stage
    x0 = [x_ode(length(x_ode),1); x_ode(length(x_ode),2)]; t0 = t_vec(end);
end

% Plot position and velocity
figure ; plot(t_vec, x_vec(:,1)) ; grid on ; title('position') ;
figure ; plot(t_vec, x_vec(:,2)) ; grid on ; title('velocity') ;


% Function that describes flight dynamics
function dx = flight_dynamics(t, x)
    y = x(1) ;
    dy = x(2) ;  
    g = 9.81 ;
    dx = [dy ;
          -g] ;
end

% % Function that describes stance dynamics
function dx = stance_dynamics(t, x)
y = x(1) ;
dy = x(2) ;

% Constants
k = 20000; % N/m
m = 80; % kg
L0 = 1; %m 
c = 50; %Ns/m
g = 9.81 ;

% Controller 
u = @(t) 0; % N

dx = [dy ;
    ((u(t) - k*(y-L0) - c*dy - m*g) / m)] ;
end

% Event function that describes when hopper hits the ground
function [value,isterminal,direction] = ground_contact(t,x)
    y = x(1);
    value = round(y - 1,2) ; % detect when y - L0 == 0
    isterminal = 1 ; % stop integration when y - L0 == 0
    direction = -1 ; % detect zero when function is decreasing
end

% Event function that describes when hopper leaves the ground from stance
function [value,isterminal,direction] = leaves_ground(t,x)
    y = x(1);
    value = round(y - 1,2) ; % detect when y - L0 == 0
    isterminal = 1 ; % stop integration when y - L0 == 0
    direction = 1 ; % detect only +ve y to -ve y transitions
end