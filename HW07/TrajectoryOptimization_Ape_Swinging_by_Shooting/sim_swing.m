function [t_vec, x_vec] = sim_swing(x, param)

    % unraveling the input variable
    tswing1 = x(1);
    theta1_0 = x(2); dtheta1_0 = x(3);
    theta2_0 = x(4); dtheta2_0 = x(5);
    aa1 = x(6:6+param.ngrid-1); % moment at anchored hand
    aa2 = x(6+param.ngrid:6+2*param.ngrid-1); % moment at elbow

    numperinterval = 10;
    % swing with one hand anchored
    x0 = [theta1_0; dtheta1_0; theta2_0; dtheta2_0; 0];
    options = odeset('reltol',1e-9,'abstol',1e-9);

    tinterval = tswing1/(param.ngrid-1); % interval between grid-points
    % you want to integrate only from grid-point to grid-point so you do not
    % incur inaccuracies by stepping over a grid-point
    
    t_vec = 0; x_vec = x0'; % matrices in which to store the time series of the state

    for countinterval=1:param.ngrid-1
        bb1 = [aa1(countinterval) aa1(countinterval+1)];
        bb2 = [aa2(countinterval) aa2(countinterval+1)];
        tspan = linspace((countinterval-1)*tinterval,(countinterval)*tinterval,numperinterval);

        [t,x] = ode45(@doublependulumodefile,tspan,x0,options,param,tinterval,bb1,bb2,tspan);
        x0 = x(end,:);
        t_vec = [t_vec; t(2:end)]; x_vec = [x_vec; x(2:end,:)];
    end
end