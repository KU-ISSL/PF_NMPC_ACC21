function [u,uk0] = control_RW(ref,H,Q_w,dt,x0,uk_,Q_v,N,M)
nx = size(x0,1);
nu = size(uk_,1);
xbar = zeros(nx+nu,N,H);

%  Create set of control particles
for i = 1:N
    xbar(:,i,1)= mvnrnd([x0;zeros(2,1)],Q_w);
    uk0(:,i,1)= xbar(5:6,i,1);
end

%% Particle Filtering
wk(:,1)   = repmat(1/N,N,1);
for j = 2:H
    for i = 1:N
        xbar(:,i,j) = sys_vir(xbar(:,i,j-1),dt,Q_w);                                    % Pass Particles through system dynamics
        uk0(:,i,j) = xbar(5:6,i,j);                                                     % Extract Control Partilces
        wk(i,j) = wk(i,j-1)*mvnpdf((ref(:,j)-M*xbar(:,i,j))',zeros(6,1)',Q_v);          % Weight Particles
    end
    wk(:,j) = wk(:,j)./sum(wk(:,j));                                                    % Normalize Weights
    %Resampling
    N_eff = 1/sum(wk(:,j).^2);                                                          % Effective number of samples
    N_T = N/1;                                                                          % Resampling Threshold
    if j~=H
        if N_eff < N_T
            [xbar(:,:,j),wk(:,j),idx] = resample(xbar(:,:,j),wk(:,j),'systematic_resampling');
        end
    end
end
%% Particle Smoothing
wsk(:,H) = wk(:,H);
% Reweighting Procedure Using Eq. (13)                    
for l = H-1:-1:1
    for i = 1:N
        for j =1:N
            pxj(j,i) = mvnpdf(xbar(:,j,l+1),sys_vir(xbar(:,i,l),dt,Q_w),Q_w)*wk(i,l);
        end
    end
    pxjn = sum(pxj,2);
    pxj = pxj./pxjn;
    pxj = pxj.*wsk(:,l+1);
    wsk(:,l) = sum(pxj,1);
end
% Extracting the Control 
u = zeros(nu,1);
for i = 1:N
    u = u + wsk(i,1)*xbar(nx+1:end,i,1);
end
end

