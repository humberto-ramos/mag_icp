clear all 
close all
% addpath("../matlab-mag-pf/")

% -- UNCOMMENT THIS SECTION FOR TRAJECTORY -- %
% trajectory = load("trajectory.mat");
% plot3(trajectory.trajectory_matrix(1,:),trajectory.trajectory_matrix(2,:),trajectory.trajectory_matrix(3,:))
% xlim([-1 1])


%% Generate model points to represent Magnetic Map
% Red Points generate from Fake Map (visualize_map.m file)

% Specify the range for plot (size of flight space)
x_range = -1.1:0.1:1.2;
y_range = -3.5:0.1:3.5;

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
% data = trajectory.trajectory_matrix;

% Create a zeros matrix for path data
data = zeros(3,N);


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

%% -- MANUAL ROTATION -- %%
v1=0; v2=0; v3=3.14;
R1=[1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
R2=[cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
R3=[cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];

R=R3*R2*R1;

data=R*data;
%% -- MANUAL ROTATION -- %%

% data(1,:)=data(1,:)+0.5;
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
[RotMat,TransVec,dataOut,res]=icp(model,data,[], [], 1, 0.5);
disp('Residual is: ')
res


%% Reference:
%
% Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
% Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2


%% -- PLOT: Model points and data points in transformed positions

% figure(2)
% plot3(model(1,:),model(2,:),model(3,:),'r.',dataOut(1,:),dataOut(2,:),dataOut(3,:),'g.')
% hold on
% % plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
% title('Transformed data points (green) and model points (red)')

