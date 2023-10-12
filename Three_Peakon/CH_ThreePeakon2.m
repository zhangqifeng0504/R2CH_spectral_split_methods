%% Nonlinear scheme for R2CH（三峰）加入黏性
%% Scheme 4
function [u] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu)
h = (xb-xa)/M; x = xa:h:xb; x = x';
ta = 0; tau = (tb-ta)/N; t = ta:tau:tb; t = t';
%% 初边值条件
u = zeros(M+1,N+1); 
c1 = 2; c2 = 1; c3 = 0.8; a = xb; x1 = -5; x2 = -3; x3 = -1; 
phi1 = @(x) c1/cosh(a/2) * cosh(x-x1) .* ( abs(x-x1) <= a/2 )...
           + c1/cosh(a/2) * cosh(a - (x-x1)) .* ( abs(x-x1) > a/2 );
phi2 = @(x) c2/cosh(a/2) * cosh(x-x2) .* ( abs(x-x2) <= a/2 )...
           + c2/cosh(a/2) * cosh(a - (x-x2)) .* ( abs(x-x2) > a/2 );       
phi3 = @(x) c3/cosh(a/2) * cosh(x-x3) .* ( abs(x-x3) <= a/2 )...
           + c3/cosh(a/2) * cosh(a - (x-x3)) .* ( abs(x-x3) > a/2 );
u(:,1) = phi1(x(:)) + phi2(x(:)) + phi3(x(:)); 
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
tol = 1e-8; eps = 1e-8;
for n = 1:N
   V = u(1:M,n); err = 1;
   while err >= tol
      epsilon_u = 0.025 .* ( abs(D2*V) >= eps/h ) + 0 .* ( abs(D2*V) < eps/h );  % very good
      B_hat = -diag(V)*D1 - D1*diag(V) + diag(D2*V)*D1 + D1*diag(D2*V);
      L = I - tau/2*(I-D2)^(-1)*B_hat - tau/2*A*D + mu*tau/2*D*D2 - tau*h/4*diag(epsilon_u)*D2;
      R = u(1:M,n);
      V1 = L \ R;  % u_half
      err = max(max(abs(V1 - V)));
      V = V1; 
   end
   u(1:M,n+1) = 2 * V1 - u(1:M,n); u(M+1,n+1) = u(1,n+1);
end

