function result = simulate()
% simulate a Raibert hopper

global h_axes body leg
global dt time x y xd yd 
global hip_torque leg_angle body_angle leg_angled body_angled
global leg_state foot_x foot_y leg_lengthd leg_length rest_leg_length

mass = 1.0;
g = 9.81;
leg_k = 200;
leg_damping = 1.0;

body_moi = 1.0;
leg_moi = 0.1;

foot_y_new = y - rest_leg_length*cos( leg_angle );
foot_x_new = x + rest_leg_length*sin( leg_angle );

if foot_y_new > 0
 leg_state = 0; % in air
 spring_force = 0;
 foot_x = foot_x_new;
 foot_y = foot_y_new;
 leg_length = rest_leg_length;
else
 if leg_state == 0
   foot_x = foot_x_new;
   foot_y = 0;
   leg_state = 1; % we are on the ground
 end;
 leg_length = sqrt( (x - foot_x)^2 + (y - foot_y)^2 );
 leg_lengthd = ((x - foot_x)*xd + (y - foot_y)*yd)/leg_length;
 spring_force = leg_k*(rest_leg_length - leg_length) - leg_damping*leg_lengthd;
end

fx = -spring_force*sin(leg_angle);
fy = spring_force*cos(leg_angle);

time = time + dt;

% force x to be zero
% xd_new = xd + dt*fx/mass;
% x = x + dt*(xd + xd_new)/2;
% xd = xd_new;
x = 0;
xd = 0;

% force y to not change
% yd_new = yd + dt*(-mass*g + fy)/mass;
% y = y + dt*(yd + yd_new)/2;
% yd = yd_new;
yd = 0;

body_angled_new = body_angled + dt*hip_torque/body_moi;
body_angle = body_angle + dt*(body_angled + body_angled_new)/2;
body_angled = body_angled_new;

if ( leg_state == 0 ) % in air
 leg_angled_new = leg_angled - dt*hip_torque/leg_moi;
 leg_angle = leg_angle + dt*(leg_angled + leg_angled_new)/2;
 leg_angled = leg_angled_new;
else  % on ground
 leg_angle = atan2( foot_x - x, y - foot_y );
 leg_angled = 0;
end

result = y;
