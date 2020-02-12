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