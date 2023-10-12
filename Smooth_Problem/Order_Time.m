%% 时间阶
clear all;  format short
% M = 40; xa = -6; xb = 6; tb = 20; a = 0.1; A = 0; mu = 0; Omega = 0;   % M=40; N=80
% M = 38; xa = -12*pi; xb = 12*pi; tb = 2; a = 4; A = 0; mu = 0; Omega = 0; % M=38; N=8
% M = 40; xa = -8; xb = 8; tb = 1; a = 0.2; A = 0; mu = 1; Omega = 73e-6; % M=40; N=8
M = 40; xa = -8; xb = 8; tb = 1; a = 1; A = 1; mu = 1; Omega = 73e-6;  % M=40;N=8
cputime = 0; tic;
for i = 1:5
    N = 8*2^(i-1);
    [u,rho] = R2CH_NonlinearScheme2(M,N,xa,xb,tb,a,A,mu,Omega);
    [U,Rho] = R2CH_NonlinearScheme2(M,2*N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr(i) = max(max(abs(u(:,:)-U(:,1:2:end))));
    rho_MaxErr(i) = max(max(abs(rho(:,:)-Rho(:,1:2:end))));
    cputime(i+1) = cputime(i) + toc;
end
cputime
u_Ord = log2(u_MaxErr(1:end-1)./u_MaxErr(2:end))
rho_Ord = log2(rho_MaxErr(1:end-1)./rho_MaxErr(2:end))