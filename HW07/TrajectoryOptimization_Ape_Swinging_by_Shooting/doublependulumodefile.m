function dstatevar = doublependulumodefile(t,x,param,tinterval,bb1,bb2,tspan)
% DOUBLEPENDULUMODEFILE This is the ODE file for a double pendulum with
% massless links, and point-masses at the end of each link (see text book 
% for the equations)

    % Extract constants
    m1 = param.mbody; m2 = param.mhand;
    L1 = param.Larm; L2 = L1; g = param.gravg;

    % Extract States
    theta1 = x(1);
    dtheta1 = x(2);
    theta2 = x(3);
    dtheta2 = x(4);

    % Compute Torque(t)
    Tor1 = bb1(1)+(bb1(2)-bb1(1))/tinterval*(t-tspan(1));
    Tor2 = bb2(1)+(bb2(2)-bb2(1))/tinterval*(t-tspan(1));

    % Compute Dynamics
    % E = - C*dq - G + B*u
    E = zeros(2,1);
    E(1) = -L2*m2*L1*sin(theta2)*dtheta1^2-L2*g*m2*cos(theta1+theta2)+Tor2;
    E(2) = L1*L2*m2*sin(theta2)*dtheta2^2 + ...
        2*L1*L2*dtheta1*m2*sin(theta2)*dtheta2-L2*g*m2*cos(theta1+theta2) ...
        -L1*g*m1*cos(theta1)-L1*g*m2*cos(theta1)+Tor1;

    D = zeros(2,2);
    D(1,1) = m2*L2^2 + m2*L1*cos(theta2)*L2;
    D(1,2) = L2^2*m2;
    D(2,1) = L1^2*m1 + L1^2*m2 + L2^2*m2 + 2*L1*L2*m2*cos(theta2);
    D(2,2) = m2*L2^2 + L1*m2*cos(theta2)*L2;
    ddq = D\E; % D\E is equivalent to inv(D)*E, but much faster to compute.

    ddtheta1 = ddq(1);
    ddtheta2 = ddq(2);

    % Compute State derivative
    dCost = Tor1^2+Tor2^2;

    % Compute dx
    dstatevar = [dtheta1; 
                 ddtheta1; 
                 dtheta2; 
                 ddtheta2; 
                 dCost];
end