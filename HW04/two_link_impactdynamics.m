function [xplus] = two_link_impactdynamics(xminus)

xplus = [-1 0 0 0;
          -2 0 0 0;
          0 0 cos(2*xminus(1)) 0;
          0 0 cos(2*xminus(1))*(1-cos(2*xminus(1))) 0] * xminus';  

end