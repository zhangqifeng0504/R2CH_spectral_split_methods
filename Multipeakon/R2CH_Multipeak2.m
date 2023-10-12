%% Nonlinear scheme for R2CH（非光滑解）无黏性项
%% Scheme 4
function [u,rho] = R2CH_Multipeak2(M,N,xa,xb,tb,A,mu,Omega)
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
%% 非线性迭代求解
tol = 1e-12; eps = 1e-5;
for n = 1:N
   V = u(1:M,n); W = rho(1:M,n); err = 1;
   while err >= tol
      W1 = rho(1:M,n) - tau/2 * D1*diag(W)*V;  % rho_half
      B_hat = -diag(V)*D1 - D1*diag(V) + diag(D2*V)*D1 + D1*diag(D2*V);
      L = I - tau/2*(I-D2)^(-1)*B_hat - tau/2*A*D + mu*tau/2*D*D2...
          - Omega*tau*(I-D2)^(-1)*diag(W1)*D1*diag(W1); 
      R = u(1:M,n) - tau*(1-2*Omega*A)/2*(I-D2)^(-1)*diag(D1*W1)*W1;
      V1 = L \ R;  % u_half
      err = max(max(abs(V1 - V)), max(abs(W1 - W)));
      V = V1; W = W1;
   end
   u(1:M,n+1) = 2 * V1 - u(1:M,n); u(M+1,n+1) = u(1,n+1);
   rho(1:M,n+1) = 2 * W1 - rho(1:M,n); rho(M+1,n+1) = rho(1,n+1);
end

