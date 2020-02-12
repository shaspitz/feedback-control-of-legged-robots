% Function to Simulate ball bouncing
function sim_ball_bounce()
     % Define time range to simulate the system
    Tspan = [0 10] ;
    t0 = 0 ; % Initial Time
    
    % Define the events function (stop integration when ball hits ground)
    options = odeset('Events', @ground_contact) ;

    % Initial state vector of the system (ball position ; ball velocity)
    x0 = [10;
          0] ;

    % Array to store data for all bounces
    t_vec = [] ; x_vec = [] ;
    
    % Simulate 5 bounces
    for j=1:5
       
        % Continuous Dynamics
        % Simulate the system until event (ball bounce) occurs
        [t x] = ode45(@ball_dynamics, t0+Tspan, x0, options) ;
        
        % Save simulation data
        t_vec = [t_vec; t] ;
        x_vec = [x_vec; x] ;
    
        % Discrete Impact Dynamics
        x_minus = x(end,:) ;
        x_plus = impact_dynamics(x_minus) ;
        
        % Initial condition / time to restart integration after bounce
        x0 = x_plus ; t0 = t_vec(end) ;
    end
    
	% Plot ball position and velocity
    figure ; plot(t_vec, x_vec(:,1)) ; grid on ; title('position') ;
    figure ; plot(t_vec, x_vec(:,2)) ; grid on ; title('velocity') ;
end

% Function that describes the continuous-time dynamics of the ball falling
function dx = ball_dynamics(t, x)
    y = x(1) ;
    dy = x(2) ;
    
    g = 9.81 ;
    
    dx = [dy ;
          -g] ;
end

% Function that describes the discrete-time impact dynamics that models ball impact
function x_plus = impact_dynamics(x_minus)
    eta = 0.8 ;
    x_plus = [x_minus(1) ;
              -eta*x_minus(2)] ;
end

% Event Function that describes when the ball hits the ground
function [value,isterminal,direction] = ground_contact(t,x)
    y = x(1) ;
    
    value = y ; % detect when y == 0
    isterminal = 1 ; % stop integration when y == 0
    direction = -1 ; % detect only +ve y to -ve y transitions
end







