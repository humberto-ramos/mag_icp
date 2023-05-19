function plot_true_and_estimate(t,xtrue,xk_updated,str_filter_name)
subplot(311);
plot(t,xk_updated(:,1),'b-');
hold on;
plot(t,xtrue(:,1),'r-');
title('X');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on;
legend('Estimate','True')

subplot(312);
plot(t,xk_updated(:,2),'b-');
hold on;
plot(t,xtrue(:,2),'r-');
title('Y');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on;
legend('Estimate','True')

subplot(313);
plot(t,xk_updated(:,3),'b-');
hold on;
plot(t,xtrue(:,3),'r-');
title('\theta');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on;
xlabel('Time(s)')
suptitle(sprintf('%s',str_filter_name))
legend('Estimate','True')