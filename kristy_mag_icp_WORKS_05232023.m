clear all 
close all
addpath("matlab-mag-pf/")
addpath("collected_trajectories/")
global partial_update
partial_update = 1;

% -- UNCOMMENT THIS SECTION FOR TRAJECTORY FROM FILE -- %
trajectory = load("trajectory_1.mat");


%% Generate model points to represent Magnetic Map
% Red Points generate from Fake Map (visualize_map.m file)
% Specify the range for plot (size of flight space)
x_range = -1.1:0.01:1.2;
y_range = -3.5:0.01:3.5;


%% Generate data points for path traveled 
% Blue data generated from Fake Map or Path traveled in Fake Map
% Import Data from file (import at top of this file)
data = trajectory.trajectory_matrix;


%% -- ARTIFICIAL DATA-- %%

sigma = 0.05;
true_path = data;

% -- GENERATE ESTIMATE PATH AS +SIGMA FROM TRUE PATH --

% WE MAY NEED A MORE SOPHISTICATED APPROACH THAT ACCOUNTS
%FOR CORRELATIONS. MAY BE EXTRACTED FROM COVARIANCE MATRIX OR NORMAL VECTOR TO
%TRAJECTORY
estim_path = true_path;
N = numel(estim_path(1, :));
n_slopes = N-1;
    %  Perpendicular Slope
    %  Theta
    %  Unit Vector cos
    %  Unit Vector sin
path_data = zeros(4,N);


%% -- CALCULATE SLOPES FROM TRUTH -- %
for i = 1: n_slopes
    % Calculate Perpendicular Slope
    y_diff = estim_path(2, i+1) - estim_path(2, i);
    x_diff = estim_path(1, i+1) - estim_path(1, i);
    slope = y_diff / x_diff;
    path_data(1, i) = -1 / slope;
    % Calculate Theta
    theta = atan(path_data(1, i));
    path_data(2, i) = theta;
    % Calculate Unit Vector
    % u = zeros(2, 1);
    % u(1,1) = cos(theta);
    % u(2,1) = sin(theta);
    path_data(3, i) = cos(theta);
    path_data(4, i) = sin(theta);
end    
path_data(3:4, N) = path_data(3:4, n_slopes);


%% CALCULATE ESTIMATED PATH AS OFFSET OF TRUTH
for i = 1: n_slopes
    % Estimated Trajectory
    estim_path(1:2, i) = estim_path(1:2, i) - 0.5*sigma*path_data(3:4, i);
end
estim_path(1:2, n_slopes+1) = estim_path(1:2, n_slopes+1) - 0.5*sigma*path_data(3:4, n_slopes);

initial_path = estim_path;

% Get New Mag Values using new plot points
for i=1:N
    estim_path(3,i) = aMap(estim_path(1,i),estim_path(2,i)); %Overwrite mag
end


%% -- CALCULATE VALUES FOR ALTERNATE TRAJECTORIES
est_plus_half = estim_path;
est_minus_half = estim_path;
est_plus = estim_path;
est_minus = estim_path;
est_plus2 = estim_path;
est_minus2 = estim_path;

for i = 1: n_slopes
    % +/- 0.5*sigma
    est_plus_half(1:2, i) = est_plus_half(1:2, i) + 0.5*sigma*path_data(3:4, i);
    est_minus_half(1:2, i) = est_minus_half(1:2, i) - 0.5*sigma*path_data(3:4, i);

     % +/- 1*sigma
    est_plus(1:2, i) = est_plus(1:2, i) + sigma*path_data(3:4, i);
    est_minus(1:2, i) = est_minus(1:2, i) - sigma*path_data(3:4, i);

     % +/- 1.5*sigma
    est_plus2(1:2, i) = est_plus2(1:2, i) + 1.5*sigma*path_data(3:4, i);
    est_minus2(1:2, i) = est_minus2(1:2, i) - 1.5*sigma*path_data(3:4, i);
end


%% Calculate Final Values for Each Trajectory
% +/- 0.5*sigma
est_plus_half(1:2, n_slopes+1) = est_plus_half(1:2, n_slopes+1) + 0.5*sigma*path_data(3:4, n_slopes);
for i=1:N
    est_plus_half(3,i) = aMap(est_plus_half(1,i),est_plus_half(2,i)); %Overwrite mag
end
est_minus_half(1:2, n_slopes+1) = est_minus_half(1:2, n_slopes+1) - 0.5*sigma*path_data(3:4, n_slopes);
for i=1:N
    est_minus_half(3,i) = aMap(est_minus_half(1,i),est_minus_half(2,i)); %Overwrite mag
end

% +/- 1*sigma
est_plus(1:2, n_slopes+1) = est_plus(1:2, n_slopes+1) + sigma*path_data(3:4, n_slopes);
for i=1:N
    est_plus(3,i) = aMap(est_plus(1,i),est_plus(2,i)); %Overwrite mag
end
est_minus(1:2, n_slopes+1) = est_minus(1:2, n_slopes+1) - sigma*path_data(3:4, n_slopes);
for i=1:N
    est_minus(3,i) = aMap(est_minus(1,i),est_minus(2,i)); %Overwrite mag
end

% +/- 1.5*sigma
est_plus2(1:2, n_slopes+1) = est_plus2(1:2, n_slopes+1) + 1.5*sigma*path_data(3:4, n_slopes);
for i=1:N
    est_plus2(3,i) = aMap(est_plus2(1,i),est_plus2(2,i)); %Overwrite mag
end
est_minus2(1:2, n_slopes+1) = est_minus2(1:2, n_slopes+1) - 1.5*sigma*path_data(3:4, n_slopes);
for i=1:N
    est_minus2(3,i) = aMap(est_minus2(1,i),est_minus2(2,i)); %Overwrite mag
end


%% ADD ERROR TO INITIAL TRAJECTORY (STARTING POINT)
error_x = -0.1;
initial_path(1,:) = initial_path(1,:) + error_x;


%% -- GENERATE REFERENCE MAP AS ESTIMATE TRAJECTORIES --
model = [est_plus_half, est_minus_half, estim_path, est_minus, est_plus, est_minus2, est_plus2];


%% -- GAUSSIAN NOISE -- 
sigma_nt = 0;
N = numel(initial_path(1,:));

initial_path(3, :) = initial_path(3, :) + sigma_nt*randn(1,N);


%% -- RUN ICP -- %%
[RotMat,TransVec,dataOut,res]=icp(model,initial_path,[], [], 1,0.000000001,true_path);


%% Reference:
%
% Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
% Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2


%% -- PLOT: Model points and data points in transformed positions
figure(3)
title("Y_Range Trajectory Fitting")
plot3(dataOut(1,:),dataOut(2,:),dataOut(3,:),'b*',true_path(1,:),true_path(2,:),true_path(3,:),'k.', estim_path(1,:),estim_path(2,:),estim_path(3,:),'g.', initial_path(1,:),initial_path(2,:),initial_path(3,:),'r.')  
% model(1,:),model(2,:),model(3,:),'b*', estim_path(1,:),estim_path(2,:),estim_path(3,:),'g.', true_path(1,:),true_path(2,:),true_path(3,:),'k.', 

legend('ICP solution', 'Truth','Estimate (Truth - 0.05sigx)','Initial condition')


