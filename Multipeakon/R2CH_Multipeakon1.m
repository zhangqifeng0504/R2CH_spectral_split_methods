%% Scheme 3 多峰模拟(加入黏性项)
function [u,rho] = R2CH_Multipeakon1(M,N,xa,xb,tb,A,mu,Omega)
h = (xb-xa)/M; x = xa:h:xb; x = x';
ta = 0; tau = (tb-ta)/N; t = ta:tau:tb; t = t';
%% 初边值条件
u = zeros(M+1,N+1); rho = zeros(M+1,N+1);
p1 = 1; p2 = -1; x1 = -5; x2 = 5;   % 波峰波谷
phi = @(x) p1 * exp(-abs(x-x1)) + p2 * exp(-abs(x-x2));
u(:,1) = phi(x(:));
rho(:,1) = 0.5;
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
%% 迭代求解
tol = 1e-8; eps = 1e-5;
for n = 1:N
   V = u(1:M,n); W = rho(1:M,n); err = 1;
   while err >= tol
      epsilon_u = 0.1  .* ( abs(D2*V) >= eps/h ) + 0 .* ( abs(D2*V) < eps/h );
      epsilon_rho = 0.1  .* ( abs(D2*W) >= eps/h ) + 0 .* ( abs(D2*W) < eps/h );
      W1 = rho(1:M,n) - tau/2 * D1*diag(W)*V + h*tau/4*diag(epsilon_rho)*D2*W;  % rho_half
      L = I + 3*tau/4*D*diag(V) + tau/4*D*diag(D1*V)*D1 - A*tau/2*D...
          - tau/2*D*D1*diag(D1*V) + tau*mu/2*D*D2...
          - Omega*tau*(I-D2)^(-1)*diag(W1)*D1*diag(W1)- tau*h/4*diag(epsilon_u)*D2;
      R = u(1:M,n) - (1-2*Omega*A)*tau/4*D*W1.^2;
      V1 = L \ R; 
      err = max(max(abs(V1-V)),max(abs(W1-W)));
      V = V1; W = W1;
   end
   u(1:M,n+1) = 2 * V1 - u(1:M,n); u(M+1,n+1) = u(1,n+1);
   rho(1:M,n+1) = 2 * W1 - rho(1:M,n); rho(M+1,n+1) = rho(1,n+1);
end




