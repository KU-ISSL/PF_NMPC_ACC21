function [u] = control_RW_C(ref,H,Q_vir,dt,x0,uk_,R_u,N,tH,tL,uk,M)
nx = size(x0,1);
nu = size(uk_,1);

uH = [3;deg2rad(35)];
uL = [-3;deg2rad(-35)];
xbar = zeros(nx+nu,N,H);

%% Particle Filtering
xk(:,:,1) = repmat(x0,1,N);
xbar(:,:,1) = [xk;uk(:,:,1)];
wk(:,1)   = repmat(1/N,N,1);
for j = 2:H
    for i = 1:N
        xbar(:,i,j) = sys_vir_C(xbar(:,i,j-1),dt,uk(:,i,j));                    % Pass Particles through system dynamics
        wk(i,j) = wk(i,j-1)*mvnpdf((ref(:,j)-M*xbar(:,i,j))',zeros(6,1)',R_u)*constraint_sat(xbar(:,i,j),uL,uH,tL,tH,ref);  % Weight Particles with Contraints Satisfaction
    end
    wk(:,j) = wk(:,j)./sum(wk(:,j));                                            % Normalize Weight Particles              
    %Resampling
    N_eff = 1/sum(wk(:,j).^2);                                                  % Effective number of samples
    N_T = N/1;                                                                  % Resampling Threshold
    if j~=H
        if N_eff < N_T
            [xbar(:,:,j),wk(:,j),idx] = resample(xbar(:,:,j),wk(:,j),'systematic_resampling');
        end
    end
end
%% Smoothing
wsk(:,H) = wk(:,H);
% Reweighting Procedure Using Eq. (13)  
for l = H-1:-1:1
    for i = 1:N
        for j =1:N
            pxj(j,i) = mvnpdf(xbar(:,j,l+1),sys_vir(xbar(:,i,l),dt,Q_vir),Q_vir)*wk(i,l);
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

