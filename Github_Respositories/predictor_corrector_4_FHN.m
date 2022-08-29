clc; clear;

% constant  
a = 0.7; b = 0.8; c = 5; I = 0.1; 

% define function handles
fv=@(t,v,u)  c * (v - u + I - (v^3) / 3);
fu=@(t,v,u)  (v - b*u + a)/c;

%initial conditions
t(1) = 0; v(1) = 0; u(1) = 0;

h=0.1; N = 1000-1;

A(1) = v(1);

for i = 1:3
    t(i+1) = t(i) + h;
    
    k1v = fv(t(i)     ,v(i)          ,u(i)          );
    k1u = fu(t(i)     ,v(i)          ,u(i)          );
    
    k2v = fv(t(i)+h/2 , v(i)+h/2*k1v , u(i)+h/2*k1u );
    k2u = fu(t(i)+h/2 , v(i)+h/2*k1v, u(i)+h/2*k1u );
    
    k3v = fv(t(i)+h/2 , v(i)+h/2*k2v , u(i)+h/2*k2u );
    k3u = fu(t(i)+h/2 , v(i)+h/2*k2v , u(i)+h/2*k2u );
    
    k4v = fv(t(i)+h , v(i)+h*k3v , u(i)+h  *k3u );
    k4u = fu(t(i)+h , v(i)+h*k3v , u(i)+h  *k3u );
    
    v(i+1) = v(i) + h/6 * (k1v + 2*k2v + 2*k3v + k4v);
    u(i+1) = u(i) + h/6 * (k1u + 2*k2u + 2*k3u + k4u);
    
    A(i+1) = v(i+1);
   
end

for i = 4:N
    
    t(i+1) = t(i) + h;
    
    partv1 = 55*fv(t(i),v(i),u(i)) - 59*fv(t(i-1),v(i-1),u(i-1)) + 37*fv(t(i-2),v(i-2),u(i-2));
    partv2 = - 9*fv(t(i-3),v(i-3),u(i-3));  
    partu1 = 55*fu(t(i),v(i),u(i)) - 59*fu(t(i-1),v(i-1),u(i-1)) + 37*fu(t(i-2),v(i-2),u(i-2));
    partu2 = - 9*fu(t(i-3),v(i-3),u(i-3));
    v(i+1) = v(i) + h*(partv1+partv2)/24;
    u(i+1) = u(i) + h*(partu1+partu2)/24;
    
    V01 = v(i+1);
    A1(i+1) = v(i+1);
           
    partv1 = 9*fv(t(i+1),v(i+1),u(i+1))+19*fv(t(i),v(i),u(i))-5*fv(t(i-1),v(i-1),u(i-1)) + fv(t(i-2),v(i-2),u(i-2));
    partu1 = 9*fu(t(i+1),v(i+1),u(i+1))+19*fu(t(i),v(i),u(i))-5*fu(t(i-1),v(i-1),u(i-1)) + fu(t(i-2),v(i-2),u(i-2));
    v(i+1) = v(i) + h*(partv1)/24;
    u(i+1) = u(i) + h*(partu1)/24;
    
    V02 = v(i+1);
    A2(i+1) = v(i+1);
    
    B(i+1) = abs((V01 - V02)/V01); % bagil hata
 
end

for k = 1:2:N
    C1((k+1)/2) = A1(k);
    C2((k+1)/2) = abs((A1(k) - A2(k))/A1(k));
end

max(B)

fig1 = figure('Position',get(0,'Screensize'));
plot(A2,'-','Color','k','linewidth',10);
grid on; hold on;
ylabel('v, [V]')
xlabel('Number of Samples');
title({'Predictor-Corrector';['I_s=',num2str(I),'    h=',num2str(h),'    v_0=',num2str(v(1)),'    u_0=',num2str(u(1))]}); grid on;
set(gca,'Fontsize',50);
saveas(fig1, 'fhn_abm.jpg');

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