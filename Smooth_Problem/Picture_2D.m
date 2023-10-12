%% 画剖面图
clear all; clc;
M = 400; N = 6400; 
xa = -12*pi; xb = 12*pi; tb = 16; a = 4; A = 0; mu = 0; Omega = 0; % Case I
% xa = -8; xb = 8; tb = 8; a = 1; A = 1; mu = 1; Omega = 73e-6;  % CaseII
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; 
t = ta:tau:tb; t = t';
[u1,rho1] = R2CH_SplittingScheme1(M,N,xa,xb,tb,a,A,mu,Omega);
[u2,rho2] = R2CH_SplittingScheme2(M,N,xa,xb,tb,a,A,mu,Omega);
[u3,rho3] = R2CH_NonlinearScheme1(M,N,xa,xb,tb,a,A,mu,Omega);
[u4,rho4] = R2CH_NonlinearScheme2(M,N,xa,xb,tb,a,A,mu,Omega);
figure(1)
plot(x,u1(:,end),'-','Color','r','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,u2(:,end),':','Color','1.00,0.07,0.65','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,u3(:,end),'.','Color','0.72,0.27,1.00','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,u4(:,end),'--','Color','0.39,0.83,0.07','Linewidth',1.5,'MarkerSize',7.5); hold on
% axis([-8 8 -0.8 0.8])
% axis([-12*pi 12*pi -1.5 1.5])
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel({'$\textbf{u}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title({'$\textbf{t=12}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 1}','\textbf{Scheme 2}','\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',14,'FontName','Times','location','Northwest');
% set(gca,'XTick',[-8:2:8]);
grid on 
set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.1);
figure(2)
plot(x,rho1(:,end),'-','Color','r','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,rho2(:,end),':','Color','1.00,0.07,0.65','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,rho3(:,end),'.','Color','0.72,0.27,1.00','Linewidth',1.5,'MarkerSize',7.5); hold on
plot(x,rho4(:,end),'--','Color','0.39,0.83,0.07','Linewidth',1.5,'MarkerSize',7.5); hold on
% axis([-8 8 0.8 2.4])
% axis([-12*pi 12*pi 1 2.6])
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel({'$\textbf{rho}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
title({'$\textbf{t=12}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 1}','\textbf{Scheme 2}','\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',14,'FontName','Times','location','Northeast');
% set(gca,'XTick',[-8:2:8]);
grid on 
set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.1);
