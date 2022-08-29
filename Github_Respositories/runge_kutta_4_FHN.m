clc; clear;

% constant  
a = 0.7; b = 0.8; c = 5; I = 1.5; 

% define function handles
fv=@(t,v,u)  c * (v - u + I - (v^3) / 3);
fu=@(t,v,u)  (v - b*u + a)/c;

%initial conditions
t(1) = 0; v(1) = 0; u(1) = 0;

%step size 
h = 0.01; N = 10000-1;

%update loop
for i=1:N
    
    t(i+1) = t(i) + h;
    
    k1v = fv(t(i)     ,v(i)          ,u(i)          );
    k1u = fu(t(i)     ,v(i)          ,u(i)          );
    
    k2v = fv(t(i)+h/2 , v(i)+h/2*k1v , u(i)+h/2*k1u );
    k2u = fu(t(i)+h/2 , v(i)+h/2*k1v , u(i)+h/2*k1u );
    
    k3v = fv(t(i)+h/2 , v(i)+h/2*k2v , u(i)+h/2*k2u );
    k3u = fu(t(i)+h/2 , v(i)+h/2*k2v , u(i)+h/2*k2u );
    
    k4v = fv(t(i)+h   , v(i)+h  *k3v , u(i)+h  *k3u );
    k4u = fu(t(i)+h   , v(i)+h  *k3v , u(i)+h  *k3u );
    
    v(i+1) = v(i) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);
    u(i+1) = u(i) + h/6 * (k1u + 2*k2u + 2*k3u + k4u);
    
end

fig1 = figure('Position',get(0,'Screensize'));
plot(v,'-','Color','k','linewidth',10);
grid on; hold on;
ylabel('v, [V]')
xlabel('Number of Samples');
title({'4th Order Runge-Kutta';['I_s=',num2str(I),'    h=',num2str(h),'    v_0=',num2str(v(1)),'    u_0=',num2str(u(1))]}); grid on;
set(gca,'Fontsize',50);
saveas(fig1, 'fhn_rk.jpg');

% fig1 = figure('Position',get(0,'Screensize'));
% plot(v,'LineStyle','-','Marker','none','Color','k','MarkerSize',20,'linewidth',4)
% grid on; hold on;
% plot(Min_locs,Min_vlue,'Marker','*','Color','r','MarkerSize',30,'linewidth',6)
% xline(Min_locs, 'Color', 'r', 'LineWidth', 4)
% textLabel = sprintf('Min of %.4f NRMSE at Polynominal Order=%.0f', Min_vlue, Min_locs);
% text(0, 130, textLabel, 'fontSize', 50, 'Color', 'r', 'VerticalAlignment','middle')
% ylabel('NRMSE')
% xlabel('Polynomial Order');
% set(gca,'Fontsize',60);
% saveas(fig1, 'fhn.jpg');

