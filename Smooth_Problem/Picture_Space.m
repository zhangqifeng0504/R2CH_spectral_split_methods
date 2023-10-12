%% 空间阶的稳定性图
clear all; clc; tic;
% N = 100; xa = -8; xb = 8; tb = 1; a = 1; A = 1; mu = 1; Omega = 73e-6;   
N = 200; xa = -12*pi; xb = 12*pi; tb = 2; a = 4; A = 0; mu = 0; Omega = 0; % CaseI 需要储存矩阵
for i = 1:8
    M(i) = 40*i;
    
    [u1,rho1] = R2CH_SplittingScheme1(M(i),N,xa,xb,tb,a,A,mu,Omega);
    [U1,Rho1] = R2CH_SplittingScheme1(2*M(i),N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr1(i) = max(max(abs(u1(:,:)-U1(1:2:end,:))));
    rho_MaxErr1(i) = max(max(abs(rho1(:,:)-Rho1(1:2:end,:))));
    
    [u2,rho2] = R2CH_NonlinearScheme2(M(i),N,xa,xb,tb,a,A,mu,Omega);
    [U2,Rho2] = R2CH_NonlinearScheme2(2*M(i),N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr2(i) = max(max(abs(u2(:,:)-U2(1:2:end,:))));
    rho_MaxErr2(i) = max(max(abs(rho2(:,:)-Rho2(1:2:end,:))));
    
    [u3,rho3] = R2CH_NonlinearScheme1(M(i),N,xa,xb,tb,a,A,mu,Omega);
    [U3,Rho3] = R2CH_NonlinearScheme1(2*M(i),N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr3(i) = max(max(abs(u3(:,:)-U3(1:2:end,:))));
    rho_MaxErr3(i) = max(max(abs(rho3(:,:)-Rho3(1:2:end,:))));
    
    [u4,rho4] = R2CH_NonlinearScheme2(M(i),N,xa,xb,tb,a,A,mu,Omega);
    [U4,Rho4] = R2CH_NonlinearScheme2(2*M(i),N,xa,xb,tb,a,A,mu,Omega);
    u_MaxErr4(i) = max(max(abs(u4(:,:)-U4(1:2:end,:))));
    rho_MaxErr4(i) = max(max(abs(rho4(:,:)-Rho4(1:2:end,:))));
end
toc
save('u_MaxErr1.mat','u_MaxErr1');
save('rho_MaxErr1.mat','rho_MaxErr1');

save('u_MaxErr2.mat','u_MaxErr2');
save('rho_MaxErr2.mat','rho_MaxErr2');

save('u_MaxErr3.mat','u_MaxErr3');
save('rho_MaxErr3.mat','rho_MaxErr3');

save('u_MaxErr4.mat','u_MaxErr4');
save('rho_MaxErr4.mat','rho_MaxErr4');

%% 画图
figure(1)
semilogy(M,u_MaxErr1(:),'g-^','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','g','MarkerEdgeColor','g')
hold on
semilogy(M,u_MaxErr2(:),'k-.v','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','k','MarkerEdgeColor','k')
hold on
semilogy(M,u_MaxErr3(:),'b--s','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','b','MarkerEdgeColor','b')
hold on
semilogy(M,u_MaxErr4(:),'r--d','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','r','MarkerEdgeColor','r')
grid on
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$M$','interpreter','latex','FontName','Times','FontSize',20);
ylabel('\textbf{Errors of} ${u}$','interpreter','latex','FontName','Times','FontSize',20);
legend('$\textbf{Scheme1}$','$\textbf{Scheme2}$','$\textbf{Scheme3}$','$\textbf{Scheme4}$',...
    'interpreter','latex','FontSize',14,'FontName','Times','Location','best')


figure(2)
semilogy(M,rho_MaxErr1(:),'g-^','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','g','MarkerEdgeColor','g')
hold on
semilogy(M,rho_MaxErr2(:),'k-.v','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','k','MarkerEdgeColor','k')
hold on
semilogy(M,rho_MaxErr3(:),'b--s','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','b','MarkerEdgeColor','b')
hold on
semilogy(M,rho_MaxErr4(:),'r--d','Linewidth',1.5,'MarkerSize',7.5,'Markerfacecolor','r','MarkerEdgeColor','r')
grid on
set(gca,'FontName','Times New Roman','FontSize',20,'FontWeight','bold','linewidth',3);
xlabel('$M$','interpreter','latex','FontName','Times','FontSize',20);
ylabel('\textbf{Errors of} ${\rho}$','interpreter','latex','FontName','Times','FontSize',20);
legend('$\textbf{Scheme1}$','$\textbf{Scheme2}$','$\textbf{Scheme3}$','$\textbf{Scheme4}$',...
    'interpreter','latex','FontSize',14,'FontName','Times','Location','best')

