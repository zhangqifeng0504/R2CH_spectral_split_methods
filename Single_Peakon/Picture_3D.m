%% 三维图 (单峰)
clear all; clc;
M = 200; N = 4000; xa = 0; xb = 20; tb = 20; A = 0; mu = 0; Omega = 0;
% M = 200; N = 4000; xa = 0; xb = 20; tb = 20; A = 1; mu = 1; Omega = 73e-6; 
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; 
t = ta:tau:tb; t = t';
[u1,rho1] = R2CH_SinglePeakon1(M,N,xa,xb,tb,A,mu,Omega);
[u2,rho2] = R2CH_SinglePeakon2(M,N,xa,xb,tb,A,mu,Omega);
u1 = real(u1); rho1 = real(rho1); u2 = real(u2); rho2 = real(rho2);
%% 画图
figure(1)
mesh(x,t,u1'); view([40 50]); 
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
ylabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
zlabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Scheme 3-u, view(40,50)}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% contourf(x,t,u1','ShowText','on')
figure(2)
mesh(x,t,rho1');view([40 50]);
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
ylabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
zlabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Scheme 3-rho, view(40,50)}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
figure(3)
mesh(x,t,u2');view([40 50]); 
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
ylabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
zlabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Scheme 4-u, view(40,50)}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
figure(4)
mesh(x,t,rho2');view([40 50]);
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
ylabel('$\textbf{t}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times'); 
zlabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title('$\textbf{Scheme 4-rho, view(40,50)}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');