%% 守恒量
clear all;  format long
% M = 60; N = 0; xa = -6; xb = 6; tb = 0; a = 0.1; A = 0; mu = 0; Omega = 0;
M = 376; N = 5120; xa = -12*pi; xb = 12*pi; tb = 20; a = 4; A = 0; mu = 0; Omega = 0; 
% M = 80; N = 10000; xa = -8; xb = 8; tb = 10; a = 0.2; A = 0; mu = 1; Omega = 73e-6; 
% M = 80; N = 2560; xa = -8; xb = 8; tb = 20; a = 1; A = 1; mu = 1; Omega = 73e-6; 
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; t = ta:tau:tb; t = t';
%% Fourier谱微分矩阵
lambda = 2*pi/(xb-xa);
row1 = zeros(M,1); row2 = zeros(M,1); I = eye(M);
row1(2:M/2) = 1:M/2-1; row1(M/2+2:M) = -M/2+1:-1;
row2 = row1; row2(M/2+1) = M/2;
Lambda1 = 1i * lambda * diag(row1);
Lambda2 = (1i * lambda * diag(row2))^2;
F = zeros(M,M); % Fourier变换矩阵
for j = 1:M
    for k = 1:M
        F(j,k) = 1 / sqrt(M) * exp(-(j-1)*(k-1)*1i*2*pi/M);
    end
end
D1 = F'*Lambda1*F; D2 = F'*Lambda2*F; D = (I-D2)^(-1) * D1;
%% Hamiltonian
[u,rho] = R2CH_NonlinearScheme1(M,N,xa,xb,tb,a,A,mu,Omega);
u0 = u(1:M,end); r0 = rho(1:M,end);   
I1 = real( h*sum(r0) );
I2 = real( h*sum(u0 + Omega * r0.^2) ); 
% E = real( h/2 * sum( u0.^2 - (D2*u0).*u0 + (1-2*Omega*A)*r0.^2 ) );
H = real( h/2 * sum( u0.^3 + u0 .* (D1*u0).^2 - A*u0.^2 + mu*u0.*(D2*u0) + u0.*r0.^2  ) );
[I1,I2,H]