function obj_sum = walk_obj(x,param)

    % Simulate walk
    [t_vec, x_vec] = sim_walk(x, param) ;
    
    alpha = x(11:14)';
    beta = x(15:18)';
    
    u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
    beta) + v_gen(s, alpha, beta));
    
    % cumulative cost
    x_TI = x_vec(end,1);
    
    % Implement the integral cost function 
    obj_sum = 0;
    for i = 1:length(t_vec)-1
    obj_sum = obj_sum + (1 / x_TI ) * norm( u(x_vec(i,:)', alpha,...
        beta) ) * (t_vec(i+1) - t_vec(i));
    end
    
end