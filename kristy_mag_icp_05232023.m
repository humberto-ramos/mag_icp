clear all 
close all
addpath("matlab-mag-pf/")
addpath("collected_trajectories/")
global partial_update
partial_update = 1;

% -- IMPORT FOR TRAJECTORY FROM FILE -- %
trajectory = load("trajectory_2.mat");


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
estim_path= true_path;
estim_path(1, :) = true_path(1, :)+sigma;
estim_path(2, :) = true_path(2, :)+sigma;


% +1*SIGMA ESTIMATE TRAJECTORY
estim_path_plus = estim_path; %Copy all. We will overwrite x and mag.
estim_path_plus(1,:) = estim_path(1,:)+sigma; %Overwrite x
estim_path_plus(2,:) = estim_path(2,:)+sigma; %Overwrite x
%Get mag data from map (aMap) at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_plus(3,i) = aMap(estim_path_plus(1,i),estim_path_plus(2,i)); %Overwrite mag
end


% -1*SIGMA ESTIMATE TRAJECTORY
estim_path_minus = estim_path; %Copy all. We will overwrite x and mag.
estim_path_minus(1,:) = estim_path(1,:)-sigma;
estim_path_minus(2,:) = estim_path(2,:)-sigma;
%Get mag data at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_minus(3,i) = aMap(estim_path_minus(1,i),estim_path_minus(2,i)); %Overwrite mag
end


% +1.5*SIGMA ESTIMATE TRAJECTORY
estim_path_plus2 = estim_path; %Copy all. We will overwrite x and mag.
estim_path_plus2(1,:) = estim_path(1,:)+1.5*sigma; %Overwrite x
estim_path_plus2(2,:) = estim_path(2,:)+1.5*sigma; 
%Get mag data at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_plus(3,i) = aMap(estim_path_plus(1,i),estim_path_plus(2,i)); %Overwrite mag
end


% -1.5*SIGMA ESTIMATE TRAJECTORY
estim_path_minus2 = estim_path; %Copy all. We will overwrite x and mag.
estim_path_minus2(1,:) = estim_path(1,:)-1.5*sigma;
estim_path_minus2(2,:) = estim_path(2,:)-1.5*sigma;
%Get mag data at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_minus(3,i) = aMap(estim_path_minus(1,i),estim_path_minus(2,i)); %Overwrite mag
end


% +0.5*SIGMA ESTIMATE TRAJECTORY
estim_path_plus_half= estim_path; %Copy all. We will overwrite x and mag.
estim_path_plus2(1,:) = estim_path(1,:)+0.5*sigma; %Overwrite x
estim_path_plus2(2,:) = estim_path(2,:)+0.5*sigma; %Overwrite x
%Get mag data at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_plus(3,i) = aMap(estim_path_plus(1,i),estim_path_plus(2,i)); %Overwrite mag
end


% -0.5*SIGMA ESTIMATE TRAJECTORY
estim_path_minus_half = estim_path; %Copy all. We will overwrite x and mag.
estim_path_minus2(1,:) = estim_path(1,:)-0.5*sigma;
estim_path_minus2(2,:) = estim_path(2,:)-0.5*sigma;
%Get mag data at this new trajectory
for i=1:numel(estim_path(1,:))
    estim_path_minus(3,i) = aMap(estim_path_minus(1,i),estim_path_minus(2,i)); %Overwrite mag
end


%% ADD ERROR TO INITIAL TRAJECTORY (STARTING POINT)
error_x = 0.1;
% estimated_path = true_path;
initial_path = estim_path;
initial_path(1,:) = initial_path(1,:) + error_x;


%% -- GENERATE REFERENCE MAP AS ESTIMATE TRAJECTORIES --
model = [estim_path_minus, estim_path_plus, estim_path_minus2, estim_path_plus2, estim_path_plus_half, estim_path_minus_half]; % estim_path];


%% -- NOISE -- %%
for i=1:numel(initial_path(1,:))
    plus_minus = round(rand);
    if plus_minus == 0
        initial_path(3, i) = initial_path(3, i) - rand*20;
    else
        initial_path(3, i) = initial_path(3, i) + rand*10;
    end
end

fprintf("Initial path: ")
initial_path(1:3, :)


%% -- Print Residuals -- %%
[RotMat,TransVec,dataOut,res]=icp(model,initial_path,[], [], 1,[],true_path);


residual = dataOut(1:2,:)-data(1:2, :);
% disp('Final End X, Y: ')
% data(1:2, end)

%% Reference:
%
% Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
% Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2


%% -- PLOT: Model points and data points in transformed positions
figure(3)
title("Y_Range Trajectory Fitting")
plot3(dataOut(1,:),dataOut(2,:),dataOut(3,:),'b*', true_path(1,:),true_path(2,:),true_path(3,:),'k.', estim_path(1,:),estim_path(2,:),estim_path(3,:),'g.',initial_path(1,:),initial_path(2,:),initial_path(3,:),'r.')  
legend('ICP solution','Truth','Estimate (Truth - 0.05sigx)','Initial condition')

