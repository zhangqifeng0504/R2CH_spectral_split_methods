%% 画图（单峰）
clear all; clc;
M = 100; N = 2000; xa = 0; xb = 20; tb = 20; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; 
t = ta:tau:tb; t = t';
[u1,rho1] = R2CH_SinglePeakon1(M,N,xa,xb,tb,A,mu,Omega);
[u2,rho2] = R2CH_SinglePeakon2(M,N,xa,xb,tb,A,mu,Omega);
figure(1)
plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title({'$\textbf{t=1}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
hold on
figure(2)
plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title({'$\textbf{t=1}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
hold on
%  set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);