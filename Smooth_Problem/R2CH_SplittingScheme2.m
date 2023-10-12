%% Splitting Method for solving R2CH
%% 空间Fourier pseudo-spectral method
function [u,rho] = R2CH_SplittingScheme2(M,N,xa,xb,tb,a,A,mu,Omega)
% M = 10; N = 10; xa = -6; xb = 6; tb = 20; a = 0.1; A = 0; mu = 0; Omega = 0;
h = (xb-xa)/M; x = xa:h:xb; x = x';
ta = 0; tau = (tb-ta)/N; t = ta:tau:tb; t = t';
%% 初边值条件
u = zeros(M+1,N+1); rho = zeros(M+1,N+1);
u(:,1) = 0; 
phi = @(x) 1 + tanh(x+a) - tanh(x-a);
rho(:,1) = phi(x(:));
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
%% 线性问题的矩阵
H0 = ( I - Omega*tau/2 * (I-D2)^(-1)*diag(rho(1:M,1))*D1*diag(rho(1:M,1)) )^(-1)...
     * ( I + Omega*tau/2 * (I-D2)^(-1)*diag(rho(1:M,1))*D1*diag(rho(1:M,1)) );
H = kron([1 0;0 0],H0) + kron([0 0;0 1],I);
B0 = -diag(u(1:M,1))*D1 - D1*diag(u(1:M,1)) + diag(D2*u(1:M,1))*D1 + D1*diag(D2*u(1:M,1));
%% 非线性问题Crank-Nicolson线性化第一层矩阵
K0 = I - tau/2*(I-D2)^(-1)*B0 - A*tau/2*D + mu*tau/2*D*D2...
     - (1-2*Omega*A)*tau^2/4*(I-D2)^(-1)*diag(rho(1:M,1))*D1*D1*diag(rho(1:M,1));
K1 = ( kron([1 0;0 0],K0) + kron([0 0;1 0],tau/2*D1*diag(rho(1:M,1))) + kron([0 0;0 1],I) )^(-1)...
    * ( kron([1 0;0 1],I) + kron([0 1;0 0],-(1-2*Omega*A)*tau/2*(I-D2)^(-1)*diag(rho(1:M,1))*D1)  );
K = 2*K1 - kron([1 0;0 1],I);
%% 第一层u1和rho1的求解
V = H*K*H * ( kron([1;0],u(1:M,1)) + kron([0;1],rho(1:M,1)) );
u(1:M,2) = V(1:M); u(M+1,2) = u(1,2); 
rho(1:M,2) = V(M+1:2*M); rho(M+1,2) = rho(1,2);
%% 求解第2~N层
for n = 2:N
   u_hat = ( 3*u(1:M,n)-u(1:M,n-1) )/2;
   rho_hat = ( 3*rho(1:M,n)-rho(1:M,n-1) )/2;    
   B_hat = -diag(u_hat)*D1 - D1*diag(u_hat) + diag(D2*u_hat)*D1 + D1*diag(D2*u_hat);
   G0 = I - tau/2*(I-D2)^(-1)*B_hat - A*tau/2*D + mu*tau/2*D*D2...
        - (1-2*Omega*A)*tau^2/4*(I-D2)^(-1)*diag(rho_hat)*D1*D1*diag(rho_hat);
   G1 = ( kron([1 0;0 0],G0) + kron([0 0;1 0],tau/2*D1*diag(rho_hat)) + kron([0 0;0 1],I) )^(-1)...
       * ( kron([1 0;0 1],I) + kron([0 1;0 0],-(1-2*Omega*A)*tau/2*(I-D2)^(-1)*diag(rho_hat)*D1) );
   G = 2*G1 - kron([1 0;0 1],I);
   W = H*G*H * ( kron([1;0],u(1:M,n)) + kron([0;1],rho(1:M,n)) );
   u(1:M,n+1) = W(1:M); u(M+1,n+1) = u(1,n+1);
   rho(1:M,n+1) = W(M+1:end); rho(M+1,n+1) = rho(1,n+1);   
end


