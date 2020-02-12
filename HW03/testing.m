clc
close all
clear all
%% Testing testing 123

% Define time range to simulate the system
Tspan = linspace(0,10,100) ;
t0 = 0 ; % Initial Time

% Initial state vector of the system (ball position ; ball velocity)
x0 = [5;
      0] ;

% Array to store data for all bounces
t_vec = [] ; x_vec = [] ;

% Simulate 5 bounces


% Define the events functions (stop integration when stance stage reached)
options1 = odeset('Events', @ground_contact) ;

% Simulate the system in flight
[t x] = ode45(@flight_dynamics, t0+Tspan, x0, options1) ;

% Save simulation data
t_vec = [t_vec; t] ;
x_vec = [x_vec; x] ;

% Initialize xo,t for stance stage
x0 = [x(length(x),1); x(length(x),2)]; t0 = t_vec(end);

% Define the events functions (stop integration when flight stage reached)
options2 = odeset('Events', @leaves_ground) ;

% Simulate the system in stance stage
[t x] = ode45(@stance_dynamics, t0+Tspan, x0, options2) ;

% Save simulation data
t_vec = [t_vec; t] ;
x_vec = [x_vec; x] ;

% Initialize xo,t for flight stage
x0 = [x(length(x),1); x(length(x),2)]; t0 = t_vec(end);

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

% Function that describes stance dynamics
function dx = stance_dynamics(t, x)
    
% Constants
k = 20000; % N/m
m = 80; % kg
L0 = 1; %m 
c = 50; %Ns/m

% Input
u = @(t) 0;

    y = x(1) ;
    dy = x(2) ;  
    g = 9.81 ;
    dx = [dy ;
          ((u(t) - k*(y-L0) - c*dy - m*g) / m)] ;
end

% Event function that describes when hopper hits the ground
function [value,isterminal,direction] = ground_contact(t,x)
    y = x(1);
    value = round(y - 1,2) ; % detect when y - L0 == 0
    isterminal = 1 ; % stop integration when y - L0 == 0
    direction = -1 ;
end

% Event function that describes when hopper leaves the ground from stance
function [value,isterminal,direction] = leaves_ground(t,x)
    y = x(1);
    value = round(y - 1,2) ; % detect when y - L0 == 0
    isterminal = 1 ; % stop integration when y - L0 == 0
    direction = 1 ; 
end