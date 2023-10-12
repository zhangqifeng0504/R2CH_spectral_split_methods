%% 画CH三峰图 t=0,1,3,4,6,8,10,12
%% t = 0
clear all;  
M = 80; N = 1000; xa = 0; xb = 30; tb = 0; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(1)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$\textbf{x}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$\textbf{u}$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',14,'FontName','Times','location','Northeast');
title({'$\textbf{t=0}$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
% set(gca,'FontName','Times New Roman','FontSize',18,'FontWeight','bold','linewidth',3); % 坐标轴设置
% grid on 
% set(gca,'GridLineStyle','-.','GridColor','k','GridAlpha',0.1); % 网格线设置
% set(gca,'XTick',[0 15 30]);
hold on
axis([0 30 -0.5 2.5])
%% t = 1
clear all;  
M = 800; N = 2000; xa = 0; xb = 30; tb = 1; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(2)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
title({'$t=1$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t =2
clear all;  
M = 800; N = 4000; xa = 0; xb = 30; tb = 2; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(3)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
title({'$t=2$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 3
clear all;  
M = 800; N = 6000; xa = 0; xb = 30; tb = 3; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(4)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northwest');
title({'$t=3$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 4
clear all;  
M = 800; N = 4000; xa = 0; xb = 30; tb = 4; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(5)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northwest');
title({'$t=4$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 6
clear all;  
M = 800; N = 12000; xa = 0; xb = 30; tb = 6; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(6)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northwest');
title({'$t=6$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 8
clear all;  
M = 800; N = 16000; xa = 0; xb = 30; tb = 8; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(7)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northwest');
title({'$t=8$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 10
clear all;  
M = 800; N = 20000; xa = 0; xb = 30; tb = 10; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(8)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
title({'$t=10$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
axis([0 30 -0.5 2.5])
%% t = 12
clear all;  
M = 800; N = 24000; xa = 0; xb = 30; tb = 12; A = 0; mu = 0;
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; t = ta:tau:tb;
[u1] = CH_ThreePeakon1(M,N,xa,xb,tb,A,mu);
[u2] = CH_ThreePeakon2(M,N,xa,xb,tb,A,mu);
figure(9)
plot(x,u1(:,end),'-','Color','0.47,0.67,0.19','Linewidth',1.5); hold on
plot(x,u2(:,end),'--','Color','1.00,0.07,0.65','Linewidth',1.5)
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel({'$x$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
ylabel('$u$','FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
legend({'\textbf{Scheme 3}','\textbf{Scheme 4}'},'interpreter','latex','FontSize',13,'FontName','Times','location','Northeast');
title({'$t=12$'},'FontUnits','points','interpreter','latex','FontSize',20,'FontName','Times');
hold on
grid on
axis([0 30 -0.5 2.5])
set(gca,'XTick',[0:5:30]);
