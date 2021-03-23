function [x1] = sys_vir(x,dt,Q_v)
L = 1.5;              % Car Length         
nx= length(x);
% Vurtual System Eq. (7) w/o Control sampling 
beta = atan(0.5*tan(x(6)));
x1 = [x(1) + x(3)*dt*cos(x(4)+beta);...
      x(2) + x(3)*dt*sin(x(4)+beta);...
      x(3) + dt*x(5);...
      x(4) + dt*x(3)*sin(beta)/(L/2);0;0];
x1 = x1 + mvnrnd(zeros(nx,1),Q_v)';
end
