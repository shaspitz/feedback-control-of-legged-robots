function [x1] = TwoLinkPoincare(x0)

% Constrain ICs to be on the poincare section (since delta perturbs it off
% the section
x0(1) = 1/2*x0(2);

% Define time range to simulate the system
Tspan = [0 15] ;
t0 = 0 ; % Initial Time

% Define the event function that stops integration at next poincare
% intersection
options = odeset('Events', @two_link_event) ;

% Simulate the system till next poincare intersection
[t_ode, x_ode] = ode45(@two_link_dynamics, t0+Tspan, x0, options) ;

% Compute point on Poincare map after one cycle (with impact dynamics
% applied)
x1 = two_link_impactdynamics(x_ode(end,:));

end