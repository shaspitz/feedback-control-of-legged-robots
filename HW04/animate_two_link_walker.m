% Animation of two link downhill walker
% Inputs:
%  t_sol: Array of time obtained from multi-step simulation
%  x_sol: Array of states obtained from multi-step simulation
%         x_sol(j, :)' is the state at time t_sol(j)
% gamma : Ground slope used for simulation
% t_I   : End of step indices into the time array.
%         t_sol(t_I(1)) = end time of first step
%         t_sol(t_I(2)) = end time of second step
%         t_sol(t_I(end)) = end time of last step
function animate_two_link_walker(t_sol, x_sol, gamma, t_I)
    if(nargin == 3)
        t_I = length(t_sol) ; % This will do "in-place animation if t_I is not provided.
    end

    % Leg length
    L = 1.5;
    
    f = figure(10) ;
    axis manual ;

    % Position of stance foot at start of simulation
    % This gets updated at end of each step to enable continous animation.
    xst = 0;
    yst = 0;

    for j=1:length(t_I)
        if j==1
            ind_start = 1 ;
        else
            ind_start = t_I(j-1)+1 ;
        end
            for k = ind_start:t_I(j)
                % Position of hip
                xH = xst-L*sin(x_sol(k,1)-gamma);
                yH = yst+L*cos(x_sol(k,1)-gamma);

                % Position of swing foot
                xsw = xH-L*sin(x_sol(k,2)-x_sol(k,1)+gamma);
                ysw = yH-L*cos(x_sol(k,2)-x_sol(k,1)+gamma);
                
                pts = [xst yst;
                       xH  yH;
                       xsw ysw] ;
                clf ;
                line([0 10.25]-1,[0 (0-10.25)*tan(gamma)],'Color',[0 0 0],'LineWidth',2); hold on ;
                plot(pts(:,1), pts(:,2), '-o') ;
                axis([-1 9 -1 2]) ; grid on ;
                drawnow ;
                pause(0.001) ;
            end
        % Update stance leg position
        xst = xsw ; yst = ysw ;
    end
end