function [x1] = sys(x,dt,u)
L = 1.5;                    % Car Length
% System Eqs. (18)
beta = atan(0.5*tan(u(2)));
x1 = [x(1) + x(3)*dt*cos(x(4)+beta);...
      x(2) + x(3)*dt*sin(x(4)+beta);...
      x(3) + dt*u(1);...
      x(4) + dt*x(3)*sin(beta)/(L/2)];
end

