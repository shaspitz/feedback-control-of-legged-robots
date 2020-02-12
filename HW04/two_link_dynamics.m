function [dx] = two_link_dynamics(t, x)

th = x(1);
phi = x(2);
dth = x(3);
dphi = x(4);

B = 0.01;
gdivl = 1;
gam = 0.01;

D = [1+2*B*(1 - cos(phi)), -B*(1 - cos(phi));
    B*(1 - cos(phi)), -B];

C = [-B*sin(phi)*(dphi^2-2*dth*dphi);
    B*dth^2*sin(phi)];

G = [B*gdivl*(sin(th-phi-gam)-sin(th-gam))-gdivl*sin(th-gam);
    B*gdivl*sin(th-phi-gam)];

ddq = inv(D)*(-G-C);

dx = [x(3:4);
      ddq];

end