% This program performs a trajectory optimization for a two-link brachiating
% robot. The goal is to find the energy-minimal way to go from one
% hand-hold to another in a given amount of time. Clearly, the programs can
% be modified to make it look like a biped problem.

% [1] M. Srinivasan. Fifteen observations on the structure of energy 
% minimizing gaits in many simple biped models. Journal of Royal Society 
% Interface, published online, June 2010. DOI: 10.1098/rsif.2009.0544
% 
% [2] M. Gomes and A. Ruina. A five-link 2D brachiating ape model with 
% life-like zero-energy-cost motions. Journal of Theoretical Biology, 
% 237, 265-278, 2005. DOI: 10.1016/j.jtbi.2005.04.014.

clear all; close all; clc; % starting with a clean slate

global var_fncount; var_fncount = 0; % this variable keeps track of function counts, so perhaps the 

% defining the various physical parameters in the problem
param.mbody = 50;  param.mhand = 0.33*param.mbody;  
param.Larm = 1;  param.gravg = 10;

% properties of actuators
param.Torqhandmax = 50; param.Torqshoulmax = 50; % 'Tor' stands for torque here.

% parameters for optimization
param.ngrid = 6; % this should be at least 2



% unknown variables to be found by optimization: initial seed for
% optimization
tswing_0 = 1.3664; theta1_0 = -2.3072; dtheta1_0 = 0; theta2_0 = -1.6271; dtheta2_0 = 0;
Torq1list_0 = 0*ones(param.ngrid,1); Torq2list_0 = 0*ones(param.ngrid,1);
x_optim0 = [tswing_0; theta1_0; dtheta1_0; theta2_0; dtheta2_0; ...
    Torq1list_0; Torq2list_0];

% running the post-processing code just to see what the initial motion
% looks like ...
swing_postprocess(x_optim0,param);
[cineq,ceq] =swing_cons(x_optim0,param)
objcost = swing_obj(x_optim0,param)


% % using MATLAB's nonlinear optimization toolbox
Aineq = []; Bineq = []; Aeq = []; Beq = []; 
LB = [0.01; -2*pi; -3; -2*pi; -3; -param.Torqhandmax*ones(param.ngrid,1); ...
    -param.Torqshoulmax*ones(param.ngrid,1)];
UB = [2; 2*pi; 3; 2*pi; 3; param.Torqhandmax*ones(param.ngrid,1); ...
    param.Torqshoulmax*ones(param.ngrid,1)];
options = optimset('display','iter','MaxFunEvals',20000,'MaxIter',20000,'diffmaxchange',1.1*1e-5, ...
    'diffminchange',1e-5);

[x_optim,obj_optim] = fmincon(@swing_obj,x_optim0,Aineq,Bineq,Aeq,Beq,LB,UB,@swing_cons,options,param);

swing_postprocess(x_optim,param);
[cineq,ceq] = swing_cons(x_optim,param)