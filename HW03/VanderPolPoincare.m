function [x1, x, t] = VanderPolPoincare(x0)

    mu = 1;
    xdot = @(t,x) [x(2); mu*(1 - x(1)^2)*x(2) - x(1)];

    % Define the events function (stop integration after one complete cycle)
    options = odeset('Events', @cycle) ;
    
    % Define time range to simulate the system
    Tspan = linspace(0,10,10^4) ;
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