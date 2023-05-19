function plot_sigma_bounds(t,state_error,covariance,str_filter_name)

subplot(311);
plot(t,state_error(:,1),'b-');
hold on;
plot(t,3.0*sqrt(covariance(:,1,1)),'r-');
plot(t,-3.0*sqrt(covariance(:,1,1)),'r-');
set(gca,'ylim',[-2 2]);
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
title('X error');
grid on;

subplot(312);
plot(t,state_error(:,2),'b-');
hold on;
plot(t,3.0*sqrt(covariance(:,2,2)),'r-');
plot(t,-3.0*sqrt(covariance(:,2,2)),'r-');
set(gca,'ylim',[-2 2]);
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
title('Y error');
grid on;

subplot(313);
plot(t,state_error(:,3),'b-');
hold on;
plot(t,3.0*sqrt(covariance(:,3,3)),'r-');
plot(t,-3.0*sqrt(covariance(:,3,3)),'r-');
set(gca,'ylim',[-2 2]);
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
title('\theta error');
xlabel('Time(s)')
suptitle(sprintf('%s',str_filter_name))
grid on;
end