clc
close all
clear all

%% ME 292B - Homework #3 Problem 2

%% Problem 2(a)
mu = 1;
xdot = @(t,x) [x(2); mu*(1 - x(1)^2)*x(2) - x(1)];

ICs = zeros(10,2);
for i = 1:size(ICs,1)
    for j = 1:size(ICs,2)
        ICs(i,j) = 10*rand(1) - 5;
    end
end

tspan = linspace(0,50,1000);
for i = 1:size(ICs,1)
    [time, xplot] = ode45(xdot,tspan,ICs(i,:));
    plot_struct(i).x1 = xplot(:,1);
    plot_struct(i).x2 = xplot(:,2);
    plot_struct(i).time = time;
end
figure()
hold on
axis([-5,5,-5,5])
title('Phase Portrait')
xlabel('x2')
ylabel('x1')
for i = 1:10
    plot(plot_struct(i).x2, plot_struct(i).x1)
end

%% Problem 2(b)

% Part i
% VanderPolPoincare is located at end of script

% Part ii
% Chose IC of [0; 4]
x0 = [0; 4];
x1 = VanderPolPoincare(x0)

% Part iii
% Obtain sequence of x's for n = 10 (x0 and x1 already obtained)
xold = x1;
xn(1,:) = x0;
xn(2,:) = x1;
for i = 3:10
    xn(i,:) = VanderPolPoincare(xold);
end

% Part iv
figure()
plot(linspace(1,10,10), xn(:,2))
title('x2(tk) vs k')
xlabel('k')
ylabel('x2(tk)')

% Show convergence to fixed point of [0 ; 2.1733]
disp(xn(end-5:end,:))
x_fixed = xn(end,:)';

%% Problem 2(c)

% Compute A numerically using finite difference
delta = 0.01;
A = zeros(2:2);
for j = 1:2
    ej = zeros(2,1);
    ej(j,1) = 1;
    A(:,j) = ( VanderPolPoincare(x_fixed + delta*ej)...
        - VanderPolPoincare(x_fixed - delta*ej) ) / (2*delta) ;
end

% Display A
A

% Find eigenvalues of A to determine stability
eigA = eig(A)
if abs(eigA) < 1
    disp('The limit cycle of the oscillator is exponentially stable')
end


% Function for 2(b)
function [x1] = VanderPolPoincare(x0)

    mu = 1;
    xdot = @(t,x) [x(2); mu*(1 - x(1)^2)*x(2) - x(1)];

    % Define the events function (stop integration after one complete cycle)
    options = odeset('Events', @cycle) ;
    
    % Define time range to simulate the system
    Tspan = linspace(0,100,10000) ;
    t0 = 0 ; % Initial Time
    
    % Simulate system
    [t x] = ode45(xdot, t0+Tspan, x0, options) ;
    
    x1 = [x(end,1);x(end,2)];
    
    function [value,isterminal,direction] = cycle(t,x)
    value = x(1) ; % detect when x1 == 0
    isterminal = 1 ; % stop integration when y == 0
    direction = 1 ; % can only approach zero while increasing (ccw)
    end
end
