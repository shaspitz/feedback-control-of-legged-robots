function [cineq,ceq] = walk_cons(x,param)

    % unraveling the input variable
    x0 = x(1:10);
    alpha = x(11:14)';
    beta = x(15:18)';
    
    u = @(s, alpha, beta) LgLfy_gen(s, alpha, beta)^-1*(-Lf2y_gen(s, alpha, ...
        beta) + v_gen(s, alpha, beta));

    % Simulate 1 step to produce (non-decision variable) paramters that will 
    % feed into nonlinear constraints
    [t_vec, x_vec] = sim_walk(x, param);
    
    u_vec = zeros(size(x_vec, 1), 2);
    for row = 1:size(x_vec, 1)
        u_vec(row,1:2) = u(x_vec(row,:)', alpha, beta);
    end
    
    for row = 1:size(x_vec, 1)
        Fst(row,:) = Fst_gen(x_vec(row,:)', u_vec(row,:)');
    end

    % inequality constraints (nonlinear)
    cineq = [];
    
    % constraint (a) - Unilateral ground constraints
    cineq = [cineq; max(-Fst(:, 2))];
    
    % constraint (b) - Friction cone constraints
    mu_s = 0.75;
    cineq = [cineq;  max( abs(Fst(:, 1)./Fst(:, 2))) - mu_s ];
    
    % constraint (c) - Desired Average Speed Constraints
    x_TI = x_vec(end,1);
    TI = t_vec(end); 
    cineq = [cineq; -x_TI / TI + param.vd];
    
    % equality constraints
    ceq = [];
    
    % constraint (d) - Periodicity constraints 
    R = [1 0 0 0 0;
    0 1 0 0 0;
    0 0 0 1 0;
    0 0 1 0 0;
    0 0 0 0 1];
    ceq = [ceq; x0(2:5,1) - R(2:end,2:end) * x_vec(end,2:5)']; % Dont include 'x' since robot moves forward after step
    ceq = [ceq; x0(6:10,1) - R * dqPlus_gen(x_vec(end,:)')];
    
end