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
function [value,isterminal,direction] = three_link_event(t,x)

value = x(3) + x(5) - pi - pi/8; % detect when phi - 2*theta == 0 (approx)
isterminal = 1 ; % stop integration when value == 0
direction = 0 ; % detect zero when function is increasing

end