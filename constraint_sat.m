function [p] = constraint_sat(x,uL,uH,xL,xH,ref)
% This function is to calculate the probability density of the constriaints
% being satisfied at time t.

% Distance to closest point on the Track
% we first figure out if the car is on the what side of the track. This is
% done by obtaining the sum of the distance form the car to the all the
% point on one side of the track

low_points = [ref(1:2,:),xL];
high_points = [ref(1:2,:),xH];
low_diff = repmat(x(1:2),1,size(low_points,2)) - low_points;
high_diff = repmat(x(1:2),1,size(high_points,2)) - high_points;
dis_low = 0;
dis_high = 0;
for i = 1:length(low_points)
    dis_low = dis_low + norm(low_diff(:,i));
    dis_high = dis_high + norm(high_diff(:,i));
end

if dis_low < dis_high
    [d,ind] = min(vecnorm(repmat(x(1:2),1,size(xL,2))-xL));
else
    [d,ind] = min(vecnorm(repmat(x(1:2),1,size(xL,2))-xH));
end

% Softplus function to captures if the constraints are satisified,
% when constraints are satisfied the barrier function should output close
% to zero.
t = 5;
alpha = 3;
barrier = @(x) (1/t)*log(1+exp(x*alpha));
g = [x(5:6)-uH;uL-x(5:6);-d];
g_bar = barrier(g);
R_eta = 0.01*eye(5);
p = mvnpdf(g_bar,zeros(length(g_bar),1),R_eta);
end

