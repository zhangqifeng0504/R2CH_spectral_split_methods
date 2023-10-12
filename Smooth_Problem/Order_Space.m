clear all; clc;
% N = 2000; xa = -6; xb = 6; tb = 20; a = 0.1; A = 0; mu = 0; Omega = 0; 
N = 200; xa = -12*pi; xb = 12*pi; tb = 2; a = 4; A = 0; mu = 0; Omega = 0;
% N = 100; xa = -8; xb = 8; tb = 1; a = 0.2; A = 0; mu = 1; Omega = 73e-6; 
% N = 100; xa = -8; xb = 8; tb = 1; a = 1; A = 1; mu = 1; Omega = 73e-6; 
for i = 1:5
    M = 32*2^(i-1);
    [u,rho] = R2CH_NonlinearScheme1(M,N,xa,xb,tb,a,A,mu,Omega);
    [U,Rho] = R2CH_NonlinearScheme1(2*M,N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr(i) = max(max(abs(u(:,:)-U(1:2:end,:))));
    rho_MaxErr(i) = max(max(abs(rho(:,:)-Rho(1:2:end,:))));
end
u_Ord = log2(u_MaxErr(1:end-1)./u_MaxErr(2:end))
rho_Ord = log2(rho_MaxErr(1:end-1)./rho_MaxErr(2:end))