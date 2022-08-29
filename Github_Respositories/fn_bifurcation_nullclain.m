clc; clear;
% okalri cizmeye yariyo
[x1, x2] = meshgrid(-8:0.1:8,-8:0.1:8);
a = 0.7; b = 0.8; c = 5; 

% I = 0.1; 
% v_fixed_point = -1.1375;
% u_fixed_point = -0.5468;

% I = 1.5; 
% v_fixed_point = 1.0324;
% u_fixed_point = 2.1656;

I = 0.5; 

dx1 = c * (x1 - x2 + I -(x1.^3)./3);
dx2 = (x1 - b .* x2 + a)./c;

r = sqrt(dx1.^2 + dx2.^2);
quiver(x1, x2, dx1./r, dx2./r, 1/2, 'LineWidth',1); % oklarin cizimi
hold on;

plot(x1(1,:), -((x1(1,:).^3)./3 - x1(1,:) - I), 'Linewidth',2);
plot(x1(1,:), -((x1(1,:) + a)./b), 'Linewidth',2);
%axis equal;
axis equal;
set(gca, 'FontSize', 20);
%axis([-8 8 -8 8]);
axis([-4 4 -4 4]);

while true
    x0 = ginput(1);
    tspan = [0 10000];
    [t, x] = ode45(@(t, x) odefcn(x, a, b, c, I), tspan, x0);
    plot(x(:,1), x(:,2), 'LineWidth', 3);
    set(gca,'Fontsize',35)
    xlabel(' v '); ylabel(' u ');
    title(['v - u Phase Portrait I_s=',num2str(I)]); grid on;
    textLabel = sprintf('UnStable');
    text(-4, -1, textLabel, 'fontSize', 35, 'Color', 'k', 'VerticalAlignment','middle')

%     plot(v_fixed_point,u_fixed_point,'Marker','*','Color','k','MarkerSize',20,'linewidth',4)
%     textLabel = sprintf('Stable point (v,u)=(%.4f,%.4f)', v_fixed_point,u_fixed_point);
%     text(-4, -1, textLabel, 'fontSize', 35, 'Color', 'k', 'VerticalAlignment','middle')

end

function dxdt = odefcn(x, a, b, c, I)

    dxdt = zeros(2,1);
    dxdt(1) = c * (x(1) - x(2) + I -(x(1)^3)/3); 
    dxdt(2) = (x(1) - b * x(2) + a)/c;
end
% 
% fig1 = figure('Position',get(0,'Screensize'));
% plot(v,'-','Color','k','linewidth',10);
% grid on; hold on;
% ylabel('v, [V]')
% xlabel('Number of Samples');
% title({'4th Order Runge-Kutta';['I_s=',num2str(I),'    h=',num2str(h),'    v_0=',num2str(v(1)),'    u_0=',num2str(u(1))]}); grid on;
% set(gca,'Fontsize',50);
% saveas(fig1, 'fhn_rk.jpg');

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
