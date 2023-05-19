%Monte Carlo runs.
close all;
clear all;

DT = 0.1;
runs = 10;
Tf = 100.0;
N = floor(Tf/DT)+1;
n_dim = 3;
%First run them all ONCE so get the plot for the estimate and associated
%uncertainty.
    [xt_ekf,xk_ekf,state_errors_ekf,P_ekf,t,Pekf5,xekf5] = ekf_robot();
    [xt_ukf,xk_ukf,state_errors_ukf,P_ukf,t,Pukf5,xukf5] = ukf_robot();
    [xt_enkf,xk_enkf,state_errors_enkf,P_enkf,t,Penk5,xenk5] = enkf_robot();
    [xt_pf,xk_pf,state_errors_pf,P_pf,t,Ppf5,xpf5,prior_pf,posterior_pf] = pf_robot();
    %Now plot the ellipses.
    %%
    figure;
%     subplot(221)
hold on
    plot_ellipse(Pekf5,xekf5,'EKF','b');
    grid on;
%     subplot(222)
    plot_ellipse(Pukf5,xukf5,'UKF','m');
    grid on;
%     subplot(223)
    plot_ellipse(Penk5,xenk5,'EnKF','r');
    grid on;
%     subplot(224)
    plot_ellipse(Ppf5,xpf5,'PF','k');
    grid on;
    legend('EKF','EKF x','UKF','UKF x','EnKF','EnKF x','PF','PF x');
    title('(x,y) estimated positions and uncertanties after 5 seconds')
    xlabel('X')
    ylabel('Y')
    set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
    grid on
      
%% Now compare the averaged estimation error over 50 Montecarlo runs.
%For this problem, we compute the root mean square squared error over the
%50 runs.
    state_errors_ekf_total = zeros(N,1);
    state_errors_ukf_total = zeros(N,1);
    state_errors_enkf_total = zeros(N,1);
    state_errors_pf_total = zeros(N,1);
for i =1 :runs
    [xt_ekf,xk_ekf,state_errors_ekf,P_ekf] = ekf_robot();
    [xt_ukf,xk_ukf,state_errors_ukf,P_ukf] = ukf_robot();
    [xt_enkf,xk_enkf,state_errors_enkf,P_enkf] = enkf_robot();
    [xt_pf,xk_pf,state_errors_pf,P_pf,t,Ppf5,xpf5] = pf_robot();
    
    state_errors_ekf_total = state_errors_ekf_total + rms(state_errors_ekf,2);
    state_errors_ukf_total = state_errors_ukf_total + rms(state_errors_ukf,2);
    state_errors_enkf_total = state_errors_enkf_total + rms(state_errors_enkf,2);
    state_errors_pf_total = state_errors_pf_total + rms(state_errors_pf,2);
     
    %get the betas trajectories for each run.
    for j = 1: N
        Pekf=reshape(P_ekf(j,:,:),3,3);
        Pukf=reshape(P_ukf(j,:,:),3,3);
        Penkf=reshape(P_enkf(j,:,:),3,3);
        Ppf=reshape(P_pf(j,:,:),3,3);

        b_ekf(i,j) =  state_errors_ekf(j,:)*inv(Pekf)*state_errors_ekf(j,:)';
        b_ukf(i,j) =  state_errors_ukf(j,:)*inv(Pukf)*state_errors_ukf(j,:)';
        b_enkf(i,j) =  state_errors_enkf(j,:)*inv(Penkf)*state_errors_enkf(j,:)';
        b_pf(i,j) =  state_errors_pf(j,:)*inv(Ppf)*state_errors_pf(j,:)';
    end
end  
%%
figure;
hold on
plot(t,mean(b_ekf,1))
plot(t,mean(b_ukf,1))
plot(t,mean(b_enkf,1))
plot(t,mean(b_pf,1))
line([0 100],[2.36 2.36])
line([0 100],[3.72 3.72])
xlabel('Time (s)')
legend('EKF','UKF','EnKF', 'PF');
title('Average NEES over 50 runs');
ylim([0,10])
grid on;

    
%% Compute RMS errors (average over 50).
state_errors_avg_ekf = state_errors_ekf_total./runs;
state_errors_avg_ukf = state_errors_ukf_total./runs;
state_errors_avg_enkf = state_errors_enkf_total./runs;
state_errors_avg_pf = state_errors_pf_total./runs;

%%
figure;
hold on
plot(t,state_errors_avg_ekf(:,1));
plot(t,state_errors_avg_ukf(:,1));
plot(t,state_errors_avg_enkf(:,1));
plot(t,state_errors_avg_pf(:,1));
legend('EKF','UKF','EnKF', 'PF');
title('RMS average error (from 50 runs)');
xlabel('Time')
ylabel('RMS error')
ylim([0 1])
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on;
grid on;

%% Get the beta averages.
% % Mean computed across dimension one.
%  bekf = mean(mean(b_ekf,2))
%  bukf = mean(mean(b_ukf,2))
% %  benkf = mean(mean(b_enkf,1))
% %  bpf = mean(mean(b_pf,1))
%  nes_bound = chi2inv(0.99,n_dim*runs)/runs
% %  nees_mean = mean(eps_ekf);

%% Plot sigma bounds and state errors for the last run.
figure;
plot_sigma_bounds(t,state_errors_ekf,P_ekf,'EKF')
figure;
plot_sigma_bounds(t,state_errors_ukf,P_ukf,'UKF')
figure;
plot_sigma_bounds(t,state_errors_enkf,P_enkf,'EnKF')
figure;
plot_sigma_bounds(t,state_errors_pf,P_pf,'PF')
%Also plot the true and estimated states.
figure;
plot_true_and_estimate(t,xt_ekf,xk_ekf,'EKF')
figure;
plot_true_and_estimate(t,xt_ukf,xk_ukf,'UKF')
figure;
plot_true_and_estimate(t,xt_enkf,xk_enkf,'EnKF')
figure;
plot_true_and_estimate(t,xt_pf,xk_pf,'PF')
%% Plot the estimated and true trajectory.
figure;
plot(xt_ekf(:,1),xt_ekf(:,2));
hold on;
plot(xk_ekf(:,1),xk_ekf(:,2));
plot(xk_ukf(:,1),xk_ukf(:,2));
plot(xk_enkf(:,1),xk_enkf(:,2));
plot(xk_pf(:,1),xk_pf(:,2));
legend('Truth','EKF','UKF','EnKF', 'PF');
title('Estimated and true trajectory (last run)');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on
grid on;
%% Plot prior before right before the first measurement.
%Get the covariances in correct shape.
%Prior covariances.
Pekf_prior=reshape(P_ekf(4,:,:),3,3);
Pukf_prior=reshape(P_ukf(4,:,:),3,3);
Penkf_prior=reshape(P_enkf(4,:,:),3,3);
Ppf_prior=reshape(P_pf(4,:,:),3,3);
%Posterior covariances.
Pekf_post=reshape(P_ekf(5,:,:),3,3);
Pukf_post=reshape(P_ukf(5,:,:),3,3);
Penkf_post=reshape(P_enkf(5,:,:),3,3);
Ppf_post=reshape(P_pf(5,:,:),3,3);

%Get the priori and posterior means.
xekf_prior = xk_ekf(4,:);
xekf_post = xk_ekf(5,:);

xukf_prior = xk_ukf(4,:);
xukf_post = xk_ukf(5,:);

xenkf_prior = xk_enkf(4,:);
xenkf_post = xk_enkf(5,:);

figure;
hold on;
subplot(221)
%Particles from pf.
hold on;
plot(xk_pf(5,1),xk_pf(5,2),'+','LineWidth',20)
plot(prior_pf(:,1),prior_pf(:,2),'.');
xlabel('X')
ylabel('Y')
title('Particle filter distribution');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(222)
plot_ellipse(Pekf_prior,xekf_prior,'EKF','b');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(223)
plot_ellipse(Pukf_prior,xukf_prior,'UKF','m');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(224)
plot_ellipse(Penkf_prior,xenkf_prior,'EnKF','r');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
% legend('Truth','EKF','UKF','EnKF', 'PF');
suptitle('Prior distributions XY')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on
xlabel('X')
ylabel('Y')
grid on;

%% Plot posterior distributions
figure;
hold on;
subplot(221)
%Particles from pf.
hold on;
plot(posterior_pf(:,1),posterior_pf(:,2),'.');
plot(xk_pf(6,1),xk_pf(6,2),'+','LineWidth',20)
xlabel('X')
ylabel('Y')
title('Particle filter distribution');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(222)
plot_ellipse(Pekf_post,xekf_post,'EKF','b');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(223)
plot_ellipse(Pukf_post,xukf_post,'UKF','m');
xlabel('X')
ylabel('Y')
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)

grid on;
subplot(224)
plot_ellipse(Penkf_post,xenkf_post,'EnKF','r');
xlabel('X')
ylabel('Y')
grid on;
% legend('Truth','EKF','UKF','EnKF', 'PF');
set(gca, 'FontSize', 15); set(gca,'LineWidth',3)
grid on
suptitle('Posterior distributions XY')
grid on;

%% Measure running times.
% [tekf,tukf,tenkf,tpf]=measure_times(runs)



