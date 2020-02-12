function swing_postprocess(x,param)
% post processing and plotting various features of the motion given some 
% initial conditions and joint torques as piecewise linear functions of 
% time 

    global var_fncount; close all;

    % unraveling the input variable
    tswing1 = x(1)
    theta1_0 = x(2) 
    dtheta1_0 = x(3)
    theta2_0 = x(4) 
    dtheta2_0 = x(5)
    aa1 = x(6:6+param.ngrid-1) % moment at anchored hand
    aa2 = x(6+param.ngrid:6+2*param.ngrid-1) % moment at elbow

    numperinterval = 30;
    % swing with one hand anchored
    statevar0 = [theta1_0; dtheta1_0; theta2_0; dtheta2_0; 0];
    options = odeset('reltol',1e-9,'abstol',1e-9);

    tinterval = tswing1/(param.ngrid-1); % interval between grid-points
    % you want to integrate only from grid-point to grid-point so you do not
    % incur inaccuracies by stepping over a grid-point

    tstore = 0; statevarstore = [statevar0']; % matrices in which to store the time series of the state

    for countinterval=1:param.ngrid-1
        bb1 = [aa1(countinterval) aa1(countinterval+1)];
        bb2 = [aa2(countinterval) aa2(countinterval+1)];
        tspan = linspace((countinterval-1)*tinterval,(countinterval)*tinterval,numperinterval);
        [tlist,statevarlist] = ode45(@doublependulumodefile,tspan,statevar0,options,param,tinterval,bb1,bb2,tspan);
        statevar0 = statevarlist(end,:);
        tstore = [tstore; tlist(2:end)]; statevarstore = [statevarstore; statevarlist(2:end,:)];
    end

    theta1store = statevarstore(:,1); dtheta1store = statevarstore(:,2);
    theta2store = statevarstore(:,3); dtheta2store = statevarstore(:,4);
    Coststore = statevarstore(:,5); % cumulative cost

    % plotting various state variables curves
    figure(1); 
    subplot(511); plot(tstore,theta1store); xlabel('t'); ylabel('theta1'); title('Simulation results');
    subplot(512); plot(tstore,dtheta1store); xlabel('t'); ylabel('theta1dot');
    subplot(513); plot(tstore,theta2store); xlabel('t'); ylabel('theta2');
    subplot(514); plot(tstore,dtheta2store); xlabel('t'); ylabel('theta2dot');
    subplot(515); plot(tstore,Coststore); xlabel('t'); ylabel('Cumulative Cost');

    % plotting the torque as a function of time
    figure(2);


    % some basic animation of the swing
    figure(2); L1 = param.Larm; L2 = param.Larm;
    plot(0,0,'ko'); hold on;
    hOA = plot([0 L1],[0 0],'r'); hold on;
    hAB = plot([L1 L1+L2],[0 0],'b');
    hA = plot(L1,0,'ro'); hB = plot(L1+L2,0,'bo');
    axis equal; axis([-1.25*(L1+L2) 1.25*(L1+L2) -1.25*(L1+L2) 1.25*(L1+L2)]);
    set(hA,'markerfacecolor','r');
    set(hB,'markerfacecolor','b');
    grid on ;

    Numtimelist = length(tstore);
    for i = 1:Numtimelist
        xA = L1*cos(theta1store(i));
        yA = L1*sin(theta1store(i));
        xB = xA + L2*cos(theta1store(i)+theta2store(i));
        yB = yA + L2*sin(theta1store(i)+theta2store(i));
        set(hOA,'xdata',[0 xA],'ydata',[0 yA]);
        set(hA,'xdata',xA,'ydata',yA);
        set(hAB,'xdata',[xA xB],'ydata',[yA yB]);
        set(hB,'xdata',xB,'ydata',yB);
        pause(0.02);
    end

    %%%% computing total energy -- it remains very close to 
    % a constant if the joint torques are zero, otherwise not quite.
    TE = zeros(Numtimelist,1);
    for i = 1:Numtimelist
        theta1 = theta1store(i);
        dtheta1 = dtheta1store(i);
        theta2 = theta2store(i);
        dtheta2 = dtheta2store(i);
        vA = [-L1*sin(theta1)*dtheta1; L1*cos(theta1)*dtheta1; 0];
        vB = vA + [-L2*sin(theta1+theta2)*(dtheta1+dtheta2); ...
            L2*cos(theta1+theta2)*(dtheta1+dtheta2); 0];
        yA = L1*sin(theta1);
        yB = yA + L2*sin(theta1+theta2);
        KE = 0.5*param.mbody*norm(vA)^2+0.5*param.mhand*norm(vB)^2;
        PE = param.mbody*param.gravg*yA+param.mhand*param.gravg*yB;
        TE(i) = KE+PE;
    end
    figure(3);
    plot(tstore,TE,'r'); ylabel('Total Energy');
    xlabel('t');
    change_in_energy = max(TE)-min(TE)
end