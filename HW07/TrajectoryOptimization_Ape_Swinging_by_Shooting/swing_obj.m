function objcost = swing_obj(x,param)
% post processing and plotting various features of the motion given some 
% initial conditions and joint torques as piecewise linear functions of 
% time 

    global var_fncount;

    % Simulate swing
    [t_vec, x_vec] = sim_swing(x, param) ;

    objcost = x_vec(end,5); % cumulative cost
    
    var_fncount = var_fncount+1;
    if mod(var_fncount,200)==0
        display(var_fncount);
        % displays the function count every 10 steps so there is something to
        % look at ...
        swing_postprocess(x,param)
    end

    if mod(var_fncount,600)==0
        save; % saves the current iteration details every once in a while.
    end
end