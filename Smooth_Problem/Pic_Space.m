clear all; close all; clc;
M = [32 64 128 256 512];

load('u_MaxErr1.mat');
load('rho_MaxErr1.mat');

load('u_MaxErr2.mat');
load('rho_MaxErr2.mat');

load('u_MaxErr3.mat');
load('rho_MaxErr3.mat');

load('u_MaxErr4.mat');
load('rho_MaxErr4.mat');


figure(1)
semilogy(M,u_MaxErr1(:),'g-^','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','g','MarkerEdgeColor','g')
hold on
semilogy(M,u_MaxErr2(:),'y-v','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','y','MarkerEdgeColor','y')
hold on
semilogy(M,u_MaxErr3(:),'b--s','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','b','MarkerEdgeColor','b')
hold on
semilogy(M,u_MaxErr4(:),'r--d','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','r','MarkerEdgeColor','r')
grid on
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$\textbf{M}$','interpreter','latex','FontName','Times','FontSize',20);
ylabel('\textbf{Numerical error of $\textbf{u}$}','interpreter','latex','FontName','Times','FontSize',20);
legend('$\textbf{Scheme1}$','$\textbf{Scheme2}$','$\textbf{Scheme3}$','$\textbf{Scheme4}$',...
    'interpreter','latex','FontSize',14,'FontName','Times','Location','best')

figure(2)
semilogy(M,rho_MaxErr1(:),'g-^','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','g','MarkerEdgeColor','g')
hold on
semilogy(M,rho_MaxErr2(:),'y-v','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','y','MarkerEdgeColor','y')
hold on
semilogy(M,rho_MaxErr3(:),'b--s','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','b','MarkerEdgeColor','b')
hold on
semilogy(M,rho_MaxErr4(:),'r--d','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','r','MarkerEdgeColor','r')
grid on
ylim([10^(-5) 10^(0)])
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$\textbf{M}$','interpreter','latex','FontName','Times','FontSize',20);
ylabel('\textbf{Numerical error of $\textbf{rho}$}','interpreter','latex','FontName','Times','FontSize',20);
legend('$\textbf{Scheme1}$','$\textbf{Scheme2}$','$\textbf{Scheme3}$','$\textbf{Scheme4}$',...
    'interpreter','latex','FontSize',14,'FontName','Times','Location','best')
