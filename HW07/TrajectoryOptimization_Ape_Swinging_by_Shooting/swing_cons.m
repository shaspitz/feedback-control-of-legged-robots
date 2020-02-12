function [cineq,ceq] = swing_cons(x,param)
% post processing and plotting various features of the motion given some 
% initial conditions and joint torques as piecewise linear functions of 
% time 

    global var_fncount;

    % unraveling the input variable
    theta1_0 = x(2); dtheta1_0 = x(3);
    theta2_0 = x(4); dtheta2_0 = x(5);

    % Simulate swing
    [t_vec, x_vec] = sim_swing(x, param) ;

    theta1_vec = x_vec(:,1); dtheta1_vec = x_vec(:,2);
    theta2_vec = x_vec(:,3); dtheta2_vec = x_vec(:,4);

    % final state
    theta1_end = theta1_vec(end); dtheta1_end = dtheta1_vec(end);
    theta2_end = theta2_vec(end); dtheta2_end = dtheta2_vec(end);

    % end-effector initial and final positions
    x_0 = param.Larm*cos(theta1_0)+param.Larm*cos(theta1_0+theta2_0);
    x_end = param.Larm*cos(theta1_end)+param.Larm*cos(theta1_end+theta2_end);
    y_0 = param.Larm*sin(theta1_0)+param.Larm*sin(theta1_0+theta2_0);
    y_end = param.Larm*sin(theta1_end)+param.Larm*sin(theta1_end+theta2_end);

    % inequality constraints
    cineq = [];
    cineq = [cineq; x_0; -x_end;];

    % equality constraints
    ceq = [];
    % constraint on initial and final hand position
    ceq = [ceq; y_0; y_end]; 
    ceq = [ceq; dtheta1_0; dtheta2_0; dtheta1_end; dtheta2_end]; % initial and final angular speeds are zero

end