%% 守恒量
clear all; clc; format long
% M = 34; N = 4000; xa = 0; xb = 1; tb = 8; A = 0; mu = 0; Omega = 0;  
M = 20; N = 4000; xa = 0; xb = 1; tb = 8; A = 0; mu = 0; Omega = 0.1;  
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
%% Hamiltonian画图
% [u1,rho1] = R2CH_Antipeak1(M,N,xa,xb,tb,A,mu,Omega);
% u_hat = u1(1:M,:); r_hat = rho1(1:M,:); 
[u2,rho2] = R2CH_Antipeak2(M,N,xa,xb,tb,A,mu,Omega);
u0 = u2(1:M,:); r0 = rho2(1:M,:); 
for n = 1:N+1
%     I1_hat(n) = real( h*sum(r_hat(:,n)) );
%     I2_hat(n) = real( h*sum(u_hat(:,n) + Omega * r_hat(:,n).^2) ); 
    I1(n) = real( h*sum(r0(:,n)) );
    I2(n) = real( h*sum(u0(:,n) + Omega * r0(:,n).^2) ); 
    E(n) = real( h/2 * sum( u0(:,n).^2 - (D2*u0(:,n)).*u0(:,n) + (1-2*Omega*A)*r0(:,n).^2 ));
end
% Err1_hat = abs(I1_hat-I1_hat(1)); Err2_hat = abs(I2_hat-I2_hat(1));
Err1 = abs(I1-I1(1)); Err2 = abs(I2-I2(1)); Err3 = abs(E-E(1)); 
% figure(1)
% plot(t,Err1_hat,'b-',t,Err1,'r-','linewidth',1.5);
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% figure(2)
% plot(t,Err2_hat,'b-',t,Err2,'r-','linewidth',1.5)
% xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% figure(3)
% plot(t,Err3,'r-','linewidth',1.5)
% xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 

figure(1)
plot(t,Err1,'r-','linewidth',1.5);
xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Mass}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
figure(2)
plot(t,Err2,'r-','linewidth',1.5)
xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Momentum}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
figure(3)
plot(t,Err3,'r-','linewidth',1.5)
xlabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{Error}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Energy}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');



