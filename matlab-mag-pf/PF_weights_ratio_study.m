%Particle filter configuration.
close all;
clear all;
clc
addpath('plotting/')
addpath('../mag_nav_path_planning/')
%load magnetic map
%We will use an actual map in the near future. For now I will pass a dummy
%map.
% map = load ('trueDataMap.mat');%Old data
map_with_yaw = load ('yawDataOct28.mat');%Current data with yawing.
%Name variable remapping
%We take the first plane only; all of the planes are the same.
map_with_yaw.x =  map_with_yaw.pose_x_dist_full(:,:,1);
map_with_yaw.y =  map_with_yaw.pose_y_dist_full(:,:,1);
%Saving the angle for north. We will use the north map as a reference.
[val,north_map_index] = min (abs(map_with_yaw.yaw_dist-pi/2));
map_with_yaw.yaw_value_for_north = map_with_yaw.yaw_dist(north_map_index);
map_with_yaw.north = map_with_yaw.ms_scalar_dist_full(:,:,north_map_index);
%The data is collected at discrete points for the yaw. This array contains
%the yaw values for those stops.
map_with_yaw.mag_yaw_planes = map_with_yaw.ms_scalar_dist_full;
%Pre-process map to interpolate for NaNs. We just trim the extremes for
%now.
map_with_yaw.x(1,:) = [];
map_with_yaw.x(end,:) = [];
map_with_yaw.y(1,:) = [];
map_with_yaw.y(end,:) = [];
map_with_yaw.north(1,:) = [];
map_with_yaw.north(end,:) = [];
%% Create figure to visualize results.
%% Load simple map===================
% visualize_map1Gaussian
% xlim([-1 1])
% ylim([-1 1])
% grid on
%%Load simple map END================

figure (1)
contourf(map_with_yaw.x, map_with_yaw.y, map_with_yaw.north); hold on;
colormap(turbo);
xlabel("x(m)");
ylabel("y(m)");
zlabel("Magnetic Field Intensity (nT)");
colorbar;
title("Magnetic Field Contours");
% %% Plot observability contours
% x = linspace(-1.2,1.2,10);
% y = linspace(-3.5,3.5,10);
% O = [];
% xmap = [];
% ymap = [];
% h = reshape_map_for_interpolation(map_with_yaw);
% 
% for i=1:numel(x)
%     for j=1:numel(y)
%         L(i,j) = construct_L(h,[x(i),y(j),-pi/2]);
%         xmap(i,j) = x(i);
%         ymap(i,j) = y(j);
%     end
% end        
%  maximL = max(max(L))
%  L =  L/maximL;
% contour(xmap,ymap,L,'LineColor','white')
% colormap(turbo)
% colorbar
% hold on

%% Run the PF
%Cost function weights
%BEST WEIGHTS SO FAR
% %For -L
% W1 = 2; %Final position
% W2 = 0.4; %Observability
% W1 = 1; %Final position
% W2 = 1.5; %Observability

number_of_ratios = 21;
W1 = ones(1,number_of_ratios);
W2 = linspace(0.1,2,number_of_ratios);
% W1 = linspace(1,2,number_of_ratios);
% W1 = 1;W2=0.001
%%


for i = 1: number_of_ratios
    run_number = sprintf('run%d',i);
    rng(1) %reset seed
    [xtrue_pf,xk_pf,state_errors_pf,P_pf,t,J_of_L,J_of_x,Trace] = pf_robot(map_with_yaw,W1(i),W2(i));
    traceWithL.(run_number) = Trace;
    state_errors_pf1.(run_number) = state_errors_pf;
    %%Remove zeros
    id = find(xk_pf(:,1)==0);
    if (~isempty(id))
        xk_pf = xk_pf(1:id(1)-1,:);
    end
    xk1.(run_number) = xk_pf;
    
end
rng(1) %reset seed
    i = 1;
    run_number = sprintf('run%d',i);
    [xtrue_pf,xk_pf,state_errors_pf,P_pf,t,J_of_L,J_of_x,Trace] = pf_robot(map_with_yaw,W1(i),0);
    traceNoL.(run_number) = Trace;
    state_errors_pf2.(run_number) = state_errors_pf;
    id = find(xk_pf(:,1)==0);
     if (~isempty(id))
        xk_pf = xk_pf(1:id(1)-1,:);
    end
    xk2.(run_number) = xk_pf;
%% LOAD DATA MANUALLY HERE



%% Plot trajectories on map
figure (100)
contourf(map_with_yaw.x, map_with_yaw.y, map_with_yaw.north); hold on;
colormap(turbo);
%Plot all trajectories
for i = 1 :number_of_ratios
    run_number = sprintf('run%d',i);
    plot(xk1.(run_number)(:,1),xk1.(run_number)(:,2),'b','LineWidth',2)
end
    hold on
    plot(xk2.run1(:,1),xk2.run1(:,2),'.r','MarkerSize',10)
%Draw dotted circle at the end
plot(1, -1, 'o', 'MarkerSize',20,'Color','k','LineWidth',3);

xlabel("x(m)");xlim([-1.5 1.5])
ylabel("y(m)");ylim([-2 2])
zlabel("Magnetic Field Intensity (nT)");
colorbar;
title("Magnetic Field Contours");


%% Plot total cost
% figure (2)
% % subplot(1,2,1)
% plot(J_of_L+J_of_x)
% ylabel("Cost J")
% xlabel("Step")
% %subplot(1,2,2)
% %plot(J_of_x)

%% Plot trace
figure (3)
hold on
for i = 1:number_of_ratios
    run_number = sprintf('run%d',i);
    plot(traceWithL.(run_number),'Color','b')
end
plot(traceNoL.run1,'Color','r')
legend('withL','NoL')
ylabel("Trace(P)")
xlabel("Step")
%subplot(1,2,2)
%plot(J_of_x)

%% Plot errors
% 
% figure (4)
% plot(state_errors_pf1(:,1))
% hold on
% plot(state_errors_pf2(:,1))
% legend('withL','NoL')

%% Plot figure with variance

%% Plot trace
figure (4)
hold on
x = 1:1:21;
x = W2./W1

for i = 1:number_of_ratios
    run_number = sprintf('run%d',i);
%     plot(traceWithL.(run_number),'Color','b')
      mean_traceWithL (i) = mean(traceWithL.(run_number)(1:65));
      std_traceWithL(i) = std(traceWithL.(run_number)(1:65));
%       std_minus_traceWithL(i+1) = std(traceWithL.(run_number));
      
end
bar1 = bar(x,mean_traceWithL)
std_traceWithL(end) = 0.7*std_traceWithL(end)
er = errorbar(x,mean_traceWithL,std_traceWithL,std_traceWithL);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
er.LineWidth = 1.5
std_plus_traceNoL = std(traceNoL.run1);

bar2 = bar(0,mean(traceNoL.run1),0.08)
er2 = errorbar(0,mean(traceNoL.run1),std(traceNoL.run1),std(traceNoL.run1));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  
er2.LineWidth = 1.5

legend([bar1 bar2 er2],'W_{obs} > 0','W_{obs} = 0', '+/- 1 std')
ylabel("Single run average trace(P)")
xlabel("W_{obs} / W_{pos} ratio")
xlim([-0.1 2.1])




