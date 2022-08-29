clc; clear;

% constant  
a = 0.7; b = 0.8; c = 5; I = 0.35; 
%a = 0.7; b = 0.8; c = 3; I = 0.34; 

% define function handles
fv=@(t,v,w)  c * (v - w + I - (v^3) / 3);
fw=@(t,v,w)  (v - b*w + a)/c;

%initial conditions
t(1) = 0; v(1) = -2; w(1) = 1.16388;

t(1) = 0; v(1) = 0; w(1) = 0;

%step size 
h=0.01; N = 10000;

%update loop
for i=1:N
    
    t(i+1) = t(i) + h;
    
    k1v = fv(t(i)     ,v(i)          ,w(i)          );
    k1w = fw(t(i)     ,v(i)          ,w(i)          );
    
    k2v = fv(t(i)+h/2 , v(i)+h/2*k1v , w(i)+h/2*k1w );
    k2w = fw(t(i)+h/2 , v(i)+h/2*k1v , w(i)+h/2*k1w );
    
    k3v = fv(t(i)+h/2 , v(i)+h/2*k2v , w(i)+h/2*k2w );
    k3w = fw(t(i)+h/2 , v(i)+h/2*k2v , w(i)+h/2*k2w );
    
    k4v = fv(t(i)+h   , v(i)+h  *k3v , w(i)+h  *k3w );
    k4w = fw(t(i)+h   , v(i)+h  *k3v , w(i)+h  *k3w );
    
    k5v = fv(t(i)+h   , v(i)+h * ((5/32)*k1v + (7/32)*k2v + (13/32)*k3v - (1/32)*k4v), w(i)+h * ((5/32)*k1w + (7/32)*k2w + (13/32)*k3w - (1/32)*k4w)); 
    k5w = fw(t(i)+h   , v(i)+h * ((5/32)*k1v + (7/32)*k2v + (13/32)*k3v - (1/32)*k4v), w(i)+h * ((5/32)*k1w + (7/32)*k2w + (13/32)*k3w - (1/32)*k4w)); 
        
    v(i+1) = v(i) + h/6 * (k1v + 2*k2v + 2*k3v + k5v);
    w(i+1) = w(i) + h/6 * (k1w + 2*k2w + 2*k3w + k5w);
    
end
4th Order Runge-Kutta

figure(1);
plot(v,'b','linewidth',4); grid on;
set(gca,'Fontsize',30)
xlabel(' Adým Sayýsý '); ylabel(' v ');
title({'Runge-Kutta New Version';['I_s=',num2str(I),'    h=',num2str(h),'    v_0=',num2str(v(1)),'    w_0=',num2str(w(1))]}); grid on;


fig1 = figure('Position',get(0,'Screensize'));
plot(NRMSE(1:80),'LineStyle','-','Marker','o','Color','k','MarkerSize',20,'linewidth',4)
grid on; hold on;
plot(Min_locs,Min_vlue,'Marker','*','Color','r','MarkerSize',30,'linewidth',6)
xline(Min_locs, 'Color', 'r', 'LineWidth', 4)
textLabel = sprintf('Min of %.4f NRMSE at Polynominal Order=%.0f', Min_vlue, Min_locs);
text(0, 130, textLabel, 'fontSize', 50, 'Color', 'r', 'VerticalAlignment','middle')
ylabel('NRMSE')
xlabel('Polynomial Order');
set(gca,'Fontsize',60);
saveas(fig1, 'fhn.jpg');

fig2 = figure('Position',get(0,'Screensize'));
plot(u,'o','Color','r','linewidth',10);
grid on; hold on;
plot(fhn_cal_w,'-.','Color','k','linewidth',10);
ylabel('u, [V]')
xlabel('Number of Samples');
set(gca,'Fontsize',60);
legend({'Original','Curve Fitting Polynomial Method'},'Location','southwest');
saveas(fig2, 'fhn_sp7_w.jpg');