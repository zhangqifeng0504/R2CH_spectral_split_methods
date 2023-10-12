%% 画图（峰值交互作用）
clear all; clc;
M = 100; N = 1000; 
xa = -20; xb = 20; tb = 1; A = 0; mu = 0; Omega = 73e-6;           % CaseIII
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; 
t = ta:tau:tb; t = t';
[u1,rho1] = R2CH_Multipeakon1(M,N,xa,xb,tb,A,mu,Omega);
[u2,rho2] = R2CH_Multipeakon2(M,N,xa,xb,tb,A,mu,Omega);
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
% %% t = 3
% clear all; clc;
% M = 400; N = 3000; 
% xa = -20; xb = 20; tb = 3; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
% h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
% x = xa:h:xb; x = x'; 
% t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Antipeakon1(M,N,xa,xb,tb,A,mu,Omega);
% [u2,rho2] = R2CH_Antipeakon2(M,N,xa,xb,tb,A,mu,Omega);
% figure(3)
% plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=3}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% figure(4)
% plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=3}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
% hold on
% %% t = 5
% clear all; clc;
% M = 400; N = 5000; 
% xa = -20; xb = 20; tb = 5; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
% h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
% x = xa:h:xb; x = x'; 
% t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Antipeakon1(M,N,xa,xb,tb,A,mu,Omega);
% [u2,rho2] = R2CH_Antipeakon2(M,N,xa,xb,tb,A,mu,Omega);
% figure(5)
% plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=5}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% figure(6)
% plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=5}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
% hold on
% %% t=6
% clear all; clc;
% M = 400; N = 6000; 
% xa = -20; xb = 20; tb = 6; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
% h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
% x = xa:h:xb; x = x'; 
% t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Antipeakon1(M,N,xa,xb,tb,A,mu,Omega);
% [u2,rho2] = R2CH_Antipeakon2(M,N,xa,xb,tb,A,mu,Omega);
% figure(7)
% plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=6}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% figure(8)
% plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=6}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% %% t=8
% clear all; clc;
% M = 400; N = 8000; 
% xa = -20; xb = 20; tb = 8; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
% h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
% x = xa:h:xb; x = x'; 
% t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Antipeakon1(M,N,xa,xb,tb,A,mu,Omega);
% [u2,rho2] = R2CH_Antipeakon2(M,N,xa,xb,tb,A,mu,Omega);
% figure(9)
% plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=8}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% figure(10)
% plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=8}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
% hold on
% %% t=10
% clear all; clc;
% M = 400; N = 10000; 
% xa = -20; xb = 20; tb = 10; A = 1; mu = 1; Omega = 73e-6;           % CaseIII
% h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
% x = xa:h:xb; x = x'; 
% t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Antipeakon1(M,N,xa,xb,tb,A,mu,Omega);
% [u2,rho2] = R2CH_Antipeakon2(M,N,xa,xb,tb,A,mu,Omega);
% figure(11)
% plot(x,u1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,u2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=10}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on
% figure(12)
% plot(x,rho1(:,end),'-','Color','0.00,0.45,0.74','MarkerSize',10,'Linewidth',1.5); hold on
% plot(x,rho2(:,end),'r--','MarkerSize',12,'Linewidth',1.5)
% xlabel('$\textbf{x}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% ylabel('$\textbf{rho}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% title({'$\textbf{t=10}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast'); 
% hold on