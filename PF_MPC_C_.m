clc; close all;clear;
nx = 4;                                     % Number of States
nu = 2;                                     % Number of Controls
rng (25)
%% Create Reference Trajectory
% Tracking a Sine wave
L = 33;                                     % Trajectory Length
dx = 0.6;
xt = dx:dx:L;
ref = zeros(6,length(xt));
ref(1,:) = xt;
ref(2,:) = 2*sin(0.2*ref(1,:));             % Reference Trajectory
tH(1,:) = ref(1,:);                         % Lower and Upper track
tH(2,:) = ref(2,:)+0.3;
tL(1,:) = ref(1,:);
tL(2,:) = ref(2,:)-0.3;
dt = 0.2;                                   % Sample Time
t = dt:dt:dt*length(ref(1,:));
T = length(t);
%% MPC
H = 4;                                      % Prediction Horizon
Ns = 100;                                   % Number of Particles
x0 = [-0.5,-0.5,3,pi/4]';                   % Initial Car Position
x = zeros(nx,T-H);
u = zeros(nu,T-H-1);                        % Control sequence for P-NMPC
u_c = zeros(nu,T-H-1);                      % Control sequence for CAP-NMPC
x(:,1)= x0;
x_c(:,1) = x0;
a = 1e-10;                                  % Numerical stabillity close to zero
Q_v = diag([0.01,0.01,a,a,a,a]);            % Covariances
Q_w = diag([a,a,a,a,0.8,0.4]);
M = diag([1,1,0,0,0,0]);                    % Augmentation Matrix
for i = 2:T-H
    % Compute Control
    tic
    if i == 2 
        [u(:,i-1),uk] = control_RW(ref(:,i-1:H+i-2),H,Q_w,dt,x0,u(:,i-1),Q_v,Ns,M);
        [u_c(:,i-1)] = control_RW_C(ref(:,i-1:H+i-2),H,Q_w,dt,x0,u_c(:,i-1),Q_v,Ns,tH(:,i-1:H+i-2),tL(:,i-1:H+i-2),uk,M);
    else
        [u(:,i-1),uk] = control_RW(ref(:,i-1:H+i-2),H,Q_w,dt,x(:,i-1),u(:,i-1),Q_v,Ns,M);
        [u_c(:,i-1)] = control_RW_C(ref(:,i-1:H+i-2),H,Q_w,dt,x_c(:,i-1),u_c(:,i-1),Q_v,Ns,tH(:,i-1:H+i-2),tL(:,i-1:H+i-2),uk,M);
    end
    % Apply to system
    
    x(:,i) = sys(x(:,i-1),dt,u(:,i-1));
    x_c(:,i) = sys(x_c(:,i-1),dt,u_c(:,i-1));
    toc
end

%% Plotting

figure(1)
plot(ref(1,1:T-H),ref(2,1:T-H),'b',tH(1,1:T-H),tH(2,1:T-H),'r--',tL(1,1:T-H),tL(2,1:T-H),'r--');
hold on;
plot(x0(1),x0(2),'o','Color','m')
plot(x(1,1:i),x(2,1:i),'k--','LineWidth',1.4);
plot(x_c(1,1:i),x_c(2,1:i),'Color','[0.4660 0.6740 0.1880]','LineWidth',1.4);
xlabel('$x$-position (m)','interpreter','latex');
ylabel('$y$-position (m)','interpreter','latex')
legend('Reference Path','Lower Boundary','Upper Boundary','Start Position','P-NMPC','CAP-NMPC')
axes('position',[.265 .20 .25 .25])
re = ref(1,:);
indexOfInterest = (re < 12) & (re > 6); % range of t near perturbation
plot(x(1,indexOfInterest),x(2,indexOfInterest),'k--','LineWidth',1.4)
hold on;
plot(tH(1,indexOfInterest),tH(2,indexOfInterest),'r--');
plot(tL(1,indexOfInterest),tL(2,indexOfInterest),'r--');
plot(ref(1,indexOfInterest),ref(2,indexOfInterest),'b');
plot(x_c(1,indexOfInterest),x_c(2,indexOfInterest),'Color','[0.4660 0.6740 0.1880]','LineWidth',1.4);
axis tight

figure(2)
plot(1:T-H-1, u(1,1:T-H-1),'k--','LineWidth',1.3);
hold on
plot(1:T-H-1,u_c(1,1:T-H-1),'Color','[0.4660 0.6740 0.1880]','LineWidth',1.3);
plot(1:T-H-1,-3*ones(1,T-H-1),'r');
plot(1:T-H-1,3*ones(1,T-H-1),'r');
ylim([-5,5])
xlabel('Time Epoch','interpreter','latex')
ylabel('Acclearation $(\mathrm{m}/\mathrm{s}^2)$','interpreter','latex')
legend('P-NMPC','CAP-NMPC','Bounds')

figure(3)
plot(1:T-H-1, rad2deg(u(2,1:T-H-1)),'k--','LineWidth',1.3);
hold on;
plot(1:T-H-1,rad2deg(u_c(2,1:T-H-1)),'Color','[0.4660 0.6740 0.1880]','LineWidth',1.3)
plot(1:T-H-1,-35*ones(1,T-H-1),'r');
plot(1:T-H-1,35*ones(1,T-H-1),'r');
xlabel('Time Epoch','interpreter','latex')
ylabel('Steering ($\circ$)','interpreter','latex')
legend('P-NMPC','CAP-NMPC','Bounds')
%% Error Calculation
err_f = sqrt(mean((x(1:2,1:i)-ref(1:2,1:i)).^2,2));
mean(err_f)
err_br = sqrt(mean((x_c(1:2,1:i)-ref(1:2,1:i)).^2,2));
mean(err_br)
%% Cost Value Calculation 
cost_un = 0;
for l = 1:length(u)-1
    cost_un = cost_un + (x(1:2,l)-ref(1:2,l))'*Q_v(1:2,1:2)^(-1)*(x(1:2,l)-ref(1:2,l)) + u(:,l)'*Q_w(5:6,5:6)^(-1)*u(:,l);
end
cost_un = cost_un + (x(1:2,end)-ref(1:2,end))'*Q_v(1:2,1:2)^(-1)*(x(1:2,end)-ref(1:2,end))
cost_c = 0;
for l = 1:length(u_c)-1
    cost_c = cost_c + (x_c(1:2,l)-ref(1:2,l))'*Q_v(1:2,1:2)^(-1)*(x_c(1:2,l)-ref(1:2,l)) + u_c(:,l)'*Q_w(5:6,5:6)^(-1)*u_c(:,l);
end
cost_c = cost_c + (x_c(1:2,end)-ref(1:2,end))'*Q_v(1:2,1:2)^(-1)*(x_c(1:2,end)-ref(1:2,end))