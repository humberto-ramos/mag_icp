clear all 
close all
addpath("matlab-mag-pf/")
addpath("collected_trajectories/")
global partial_update
partial_update = 1;

% -- UNCOMMENT THIS SECTION FOR TRAJECTORY -- %
trajectory = load("trajectory_2.mat");
% plot3(trajectory.trajectory_matrix(1,:),trajectory.trajectory_matrix(2,:),trajectory.trajectory_matrix(3,:))
% xlim([-1 1])


%% Generate model points to represent Magnetic Map
% Red Points generate from Fake Map (visualize_map.m file)

% Specify the range for plot (size of flight space)
% x_range = -1.1:0.01:1.2;
% y_range = -3.5:0.01:3.5;
% x_range = -1.1:0.01:1.2;
% y_range = -3.5:0.01:3.5;
% x_range = -1.1:0.01:1.2;
% y_range = -3.5:0.01:3.5;

% Create zeros matrix for final model
model=zeros(3,numel(x_range)*numel(y_range));

% Initialize 2D array to store Magnetic Reading [height(Z)] values
ms_fake = zeros(numel(y_range),numel(x_range));

% Compute intensities for the x and y ranges.
k = 0;

for j=1:length(x_range)

    for i=1:length(y_range)
        k = k + 1;
        ms_fake(i,j) = aMap(x_range(j), y_range(i));
        model(1,k) =  x_range(j);
        model(2,k) =  y_range(i);
        model(3,k) =  ms_fake(i,j);   

    end
end


%% Generate data points for path traveled 
% Blue data generated from Fake Map or Path traveled in Fake Map

% Use y_range for a y set of data, and x_range for x set
N = numel(y_range);
% N = numel(x_range);

% Import Data
data = trajectory.trajectory_matrix;

% disp('Start Matrix X, Y: ')
% data(1:2, 1:end)

% disp('Orig End X, Y: ')
% data(1:2, end)

% Create a zeros matrix for path data
% data = zeros(3,N);


%% If 3 sided figure used, this will set random values within each plane
% ----- DO NOT USE AT THIS TIME -----%
% data(1:2,1:N1)=rand(2,N1);
% data([1,3],(N1+1):(N1+N2))=rand(2,N2);
% data(2,(N1+1):(N1+N2))=1;
% data(2:3,(N1+N2+1):N)=rand(2,N3);
% data(1,(N1+N2+1):N)=1;
% ----- DO NOT USE AT THIS TIME -----%


%% Transform data points to their start positions
% Data set is equal to one y value of map data
% Starts at bottom left corner, and proceeds for N points.
% data = model(1:3,1:N);

% Line across the top of the graph
% data = model(1:3,71:N:end);

% Data set is equal to one x value of map data
% data = model(1:3, 1:71:end)

% Diagonal Line -- NOT SUCCESSFUL
% data = model(1:3, 1:72:end/2)

% Diagonal Line -- full length
% data = model(1:3, 1:72:end)

% Diagonal Line 
% data = model(1:3, 1:74:end);

% Diagonal Line -- NOT SUCCESSFUL
% data = model(1:3, 1:N:end/2)

% Diagonal Line -- through center
% data = model(1:3, 100:N+1:end-1000)

% Sets of Diagonals
% data = model(1:3, 1:10:end/2)

% v1=0.6*(2*rand-1); v2=0.6*(2*rand-1); v3=0.6*(2*rand-1);


%% -- ARTIFICIAL DATA-- %%
sigma = 0.05;
true_path = data;
%Plus sigma data. WE MAY NEED A MORE SOPHISTICATED APPROACH THAT ACCOUNTS
%FOR CORRELATIONS. MAY BE EXTRACTED FROM COVARIANCE MATRIX OR NORMAL VECTOR TO
%TRAJECTORY
true_path_plus = true_path; %Copy all. We will overwrite x and mag.
true_path_plus(1,:) = true_path(1,:)+sigma; %Overwrite x
%Get mag data at this new trajectory
for i=1:numel(true_path(1,:))
    true_path_plus(3,i) = aMap(true_path_plus(1,i),true_path_plus(2,i)); %Overwrite mag
end

true_path_minus = true_path; %Copy all. We will overwrite x and mag.
true_path_minus(1,:) = true_path(1,:)-sigma;
%Get mag data at this new trajectory
for i=1:numel(true_path(1,:))
    true_path_minus(3,i) = aMap(true_path_minus(1,i),true_path_minus(2,i)); %Overwrite mag
end



error_x = 0.1;
estimated_path = true_path;
estimated_path(1,:) = estimated_path(1,:) + error_x;

%% -- MANUAL ROTATION -- %%
v1=0; v2=0; v3=0.2;
R1=[1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
R2=[cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
R3=[cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];

R=R3*R2*R1;

estimated_path=R*estimated_path;

% % initial_path_plus(1,:) = estimated_path(1,:)+sigma;
% % initial_path_minus(1,:) = estimated_path(1,:)-sigma;

model = [true_path, true_path_minus,true_path_plus];

% data(2,:)=data(2,:)+0.1;
% data(1,:)=data(1,:)+abs(0.2*randn);
% data(2,:)=data(2,:)+abs(0.2*randn);
% data(3,:)=data(3,:)+0.2*randn;

% A plot. Model points and data points in start positions
%%
% figure(1)
% % plot3(model(1,:),model(2,:),model(3,:),'r.')
% plot3(model(1,:),model(2,:),model(3,:),'r.',data(1,:),data(2,:),data(3,:),'b.')
% hold on
% 
% xlabel("X")
% ylabel("Y")
% zlabel("Z")
% title('Original data points (blue) and model points (red)')
% grid on
% Running the ICP-algorithm. Least squares criterion

%% -- Print Residuals -- %%
[RotMat,TransVec,dataOut,res]=icp(model,estimated_path,[], [], 1,[],true_path);
% disp('Residual is: ')
% res

% disp('Final Matrix X, Y: ')
% dataOut(1:2, 1:end)

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
plot3(dataOut(1,:),dataOut(2,:),dataOut(3,:),'b*', true_path(1,:),true_path(2,:),true_path(3,:),'k.', estimated_path(1,:),estimated_path(2,:),estimated_path(3,:),'r.')  
legend("ICP solution",'Truth','Initial condition')
% model(1,:),model(2,:),model(3,:),'g.', 
% dataOut(1,:),dataOut(2,:),dataOut(3,:),'b.',
% figure(2)
% plot3(model(1,:),model(2,:),model(3,:),'r.',dataOut(1,:),dataOut(2,:),dataOut(3,:),'g.')
% hold on
% % plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
% title('Transformed data points (green) and model points (red)')

