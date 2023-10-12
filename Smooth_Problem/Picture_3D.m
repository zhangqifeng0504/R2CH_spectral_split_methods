%% 三维图
clear all; clc;
M = 32; N = 800; xa = -12*pi; xb = 12*pi; tb = 50; a = 4; A = 0; mu = 1; Omega = 73e-06;
% M = 100; N = 6400; xa = -12*pi; xb = 12*pi; tb = 50; a = 4; A = 0; mu = 1; Omega = 73e-06; 
h = (xb-xa)/M; ta = 0; tau = (tb-ta)/N;
x = xa:h:xb; x = x'; 
t = ta:tau:tb; t = t';
% [u1,rho1] = R2CH_Scheme1(M,N,xa,xb,tb,a,A,mu,Omega);
% [u2,rho2] = R2CH_Scheme2(M,N,xa,xb,tb,a,A,mu,Omega);
[u3,rho3] = R2CH_NonlinearScheme1(M,N,xa,xb,tb,a,A,mu,Omega);
[u4,rho4] = R2CH_NonlinearScheme2(M,N,xa,xb,tb,a,A,mu,Omega);
% u1 = real(u1); rho1 = real(rho1);u2 = real(u2); rho2 = real(rho2); 
u3 = real(u3); rho3 = real(rho3); u4 = real(u4); rho4 = real(rho4);
%% 画图
% figure(1)
% mesh(x,t,u1');view([45 70]); % view([0 90])
% xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% zlabel('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% title('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% % axis([-12*pi 12*pi 0 50])
% figure(2)
% mesh(x,t,rho1');view([45 70]);
% xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% zlabel('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% title('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% % axis([-12*pi 12*pi 0 50])
% figure(3)
% mesh(x,t,u2');view([45 70]); % view([0 90])
% xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% zlabel('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% title('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% axis([-12*pi 12*pi 0 50])
% figure(4)
% mesh(x,t,rho2');view([45 70]);
% xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
% zlabel('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% title('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% % axis([-12*pi 12*pi 0 50])
figure(5)
mesh(x,t,u3');view([45 70]); % view([0 90])
xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
zlabel('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
title('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% axis([-12*pi 12*pi 0 50])
figure(6)
mesh(x,t,rho3');view([45 70]);
xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
zlabel('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
title('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% axis([-12*pi 12*pi 0 50])
figure(7)
mesh(x,t,u4');view([45 70]); % view([0 90])
xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
zlabel('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
title('$u(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% axis([-12*pi 12*pi 0 50])
figure(8)
mesh(x,t,rho4');view([45 70]);
xlabel('$x$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
ylabel('$t$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times'); 
zlabel('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
title('$\rho(x,t)$','FontUnits','points','interpreter','latex','FontSize',16,'FontName','Times');
% axis([-12*pi 12*pi 0 50])