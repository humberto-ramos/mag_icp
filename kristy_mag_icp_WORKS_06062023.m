clear all 
close all
addpath("matlab-mag-pf/")
addpath("collected_trajectories/")
global partial_update
partial_update = 1;

% -- UNCOMMENT THIS SECTION FOR TRAJECTORY FROM FILE -- %
traj_name = "trajectory_2";
trajectory = load(traj_name+".mat");

%% -- NOT CURRENTLY USED --
% NOT USED FOR 06/2023 WORKING VERSION
% plot3(trajectory.trajectory_matrix(1,:),trajectory.trajectory_matrix(2,:),trajectory.trajectory_matrix(3,:))
% xlim([-1 1])


%% Generate model points to represent Magnetic Map
% Red Points generate from Fake Map (visualize_map.m file)
% Specify the range for plot (size of flight space)
x_range = -1.1:0.05:1.2;
y_range = -3.5:0.05:3.5;

% %% -- NOT CURRENTLY USED --
% NOT USED FOR 06/2023 WORKING VERSION
% x_range = -1.1:0.01:1.2;
% y_range = -3.5:0.01:3.5;
% x_range = -1.1:0.01:1.2;
% y_range = -3.5:0.01:3.5;
% 
% % Create zeros matrix for final model
% model=zeros(3,numel(x_range)*numel(y_range));
% 
% % Initialize 2D array to store Magnetic Reading [height(Z)] values
% ms_fake = zeros(numel(y_range),numel(x_range));
% 
% % Compute intensities for the x and y ranges.
% k = 0;
% 
% for j=1:length(x_range)
% 
%     for i=1:length(y_range)z
%         k = k + 1;
%         ms_fake(i,j) = aMap(x_range(j), y_range(i));
%         model(1,k) =  x_range(j);
%         model(2,k) =  y_range(i);
%         model(3,k) =  ms_fake(i,j);   
% 
%     end
% end


%% Generate data points for path traveled 
% Blue data generated from Fake Map or Path traveled in Fake Map
% Import Data from file (import at top of this file)
data = trajectory.trajectory_matrix;


%% -- NOT CURRENTLY USED --
% NOT USED FOR 06/2023 WORKING VERSION
% Use y_range for a y set of data, and x_range for x set
% N = numel(y_range);
% N = numel(x_range); % Do not use

% disp('Start Matrix X, Y: ')
% data(1:2, 1:end)

% disp('Orig End X, Y: ')
% data(1:2, end)

% Create a zeros matrix for path data
% data = zeros(3,N);


%% If 3 sided figure used, this will set random values within each plane
% NOT USED FOR 06/2023 WORKING VERSION
% ----- DO NOT USE AT THIS TIME -----%
% data(1:2,1:N1)=rand(2,N1);
% data([1,3],(N1+1):(N1+N2))=rand(2,N2);
% data(2,(N1+1):(N1+N2))=1;
% data(2:3,(N1+N2+1):N)=rand(2,N3);
% data(1,(N1+N2+1):N)=1;
% ----- DO NOT USE AT THIS TIME -----%


%% -- NOT CURRENTLY USED -- 
% NOT USED FOR 06/2023 WORKING VERSION
% Transform data points to their start positions
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
initial_path = true_path;


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


%% We must update our true path to have map values.
% Here we update the true_path to have clean mag data
for i=1:N
    true_path(3,i) = aMap(true_path(1,i),true_path(2,i)); %Overwrite mag
end


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


%% CALCULATE ESTIMATED PATH AS OFFSET OF INITIAL PATH
% NOT USED FOR 06/2023 WORKING VERSION
% After this section, the initial path will have the collected values (with
% noise), with the "filter estimate" x/y values. 
% for i = 1: n_slopes
%     % Estimated Trajectory
%     initial_path(1:2, i) = initial_path(1:2, i) - 0.5*sigma*path_data(3:4, i);
% end
% initial_path(1:2, n_slopes+1) = initial_path(1:2, n_slopes+1) - 0.5*sigma*path_data(3:4, n_slopes);


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


%% -- OLD METHOD FOR MOVING POINTS --
% NOT USED FOR 06/2023 WORKING VERSION
% estim_path(1, :) = true_path(1, :)+sigma;
% % estim_path(2, :) = true_path(2, :)+sigma;
% 
% 
% % +1*SIGMA ESTIMATE TRAJECTORY
% estim_path_plus = estim_path; %Copy all. We will overwrite x and mag.
% estim_path_plus(1,:) = estim_path_plus(1,:)+sigma; %Overwrite x
% % estim_path_plus(2,:) = estim_path_plus(2,:)+sigma; %Overwrite x
% %Get mag data from map (aMap) at this new trajectory
% for i=1:numel(estim_path_plus(1,:))
%     estim_path_plus(3,i) = aMap(estim_path_plus(1,i),estim_path_plus(2,i)); %Overwrite mag
% end
% 
% 
% % -1*SIGMA ESTIMATE TRAJECTORY
% estim_path_minus = estim_path; %Copy all. We will overwrite x and mag.
% estim_path_minus(1,:) = estim_path_minus(1,:)-sigma;
% % estim_path_minus(2,:) = estim_path_minus(2,:)-sigma;
% %Get mag data at this new trajectory
% for i=1:numel(estim_path_minus(1,:))
%     estim_path_minus(3,i) = aMap(estim_path_minus(1,i),estim_path_minus(2,i)); %Overwrite mag
% end
% 
% % +1.5*SIGMA ESTIMATE TRAJECTORY
% estim_path_plus2 = estim_path; %Copy all. We will overwrite x and mag.
% estim_path_plus2(1,:) = estim_path_plus2(1,:)+1.5*sigma; %Overwrite x
% % estim_path_plus2(2,:) = estim_path_plus2(2,:)+1.5*sigma; 
% %Get mag data at this new trajectory
% for i=1:numel(estim_path_plus2(1,:))
%     estim_path_plus2(3,i) = aMap(estim_path_plus2(1,i),estim_path_plus2(2,i)); %Overwrite mag
% end
% 
% % -1.5*SIGMA ESTIMATE TRAJECTORY
% estim_path_minus2 = estim_path; %Copy all. We will overwrite x and mag.
% estim_path_minus2(1,:) = estim_path_minus2(1,:)-1.5*sigma;
% % estim_path_minus2(2,:) = estim_path_minus2(2,:)-1.5*sigma;
% %Get mag data at this new trajectory
% for i=1:numel(estim_path_minus2(1,:))
%     estim_path_minus2(3,i) = aMap(estim_path_minus2(1,i),estim_path_minus2(2,i)); %Overwrite mag
% end
% 
% % +0.5*SIGMA ESTIMATE TRAJECTORY
% estim_path_plus_half= estim_path; %Copy all. We will overwrite x and mag.
% estim_path_plus_half(1,:) = estim_path_plus_half(1,:)+0.5*sigma; %Overwrite x
% % estim_path_plus_half(2,:) = estim_path_plus_half(2,:)+0.5*sigma; %Overwrite x
% %Get mag data at this new trajectory
% for i=1:numel(estim_path_plus_half(1,:))
%     estim_path_plus_half(3,i) = aMap(estim_path_plus_half(1,i),estim_path_plus_half(2,i)); %Overwrite mag
% end
% 
% % -0.5*SIGMA ESTIMATE TRAJECTORY
% estim_path_minus_half = estim_path; %Copy all. We will overwrite x and mag.
% estim_path_minus_half(1,:) = estim_path_minus_half(1,:)-0.5*sigma;
% % estim_path_minus_half(2,:) = estim_path_minus_half(2,:)-0.5*sigma;
% %Get mag data at this new trajectory
% for i=1:numel(estim_path_minus_half(1,:))
%     estim_path_minus_half(3,i) = aMap(estim_path_minus_half(1,i),estim_path_minus_half(2,i)); %Overwrite mag
% end






%% -----NOT CURRENTLY USED---------------------
% NOT USED FOR 06/2023 WORKING VERSION
% Actual path +- 0.05 sigma
% true_path_plus = true_path; %Copy all. We will overwrite x and mag.
% true_path_plus(1,:) = true_path(1,:)+sigma; %Overwrite x
% %Get mag data at this new trajectory
% for i=1:numel(true_path(1,:))
%     true_path_plus(3,i) = aMap(true_path_plus(1,i),true_path_plus(2,i)); %Overwrite mag
% end
% 
% true_path_minus = true_path; %Copy all. We will overwrite x and mag.
% true_path_minus(1,:) = true_path(1,:)-sigma;
% %Get mag data at this new trajectory
% for i=1:numel(true_path(1,:))
%     true_path_minus(3,i) = aMap(true_path_minus(1,i),true_path_minus(2,i)); %Overwrite mag
% end
% ---------------------------

%% ADD ERROR TO INITIAL TRAJECTORY (STARTING POINT)
% NOT USED FOR 06/2023 WORKING VERSION
% error_x = 0.1;
% % estimated_path = true_path;
% % initial_path = estim_path;
% initial_path(1,:) = initial_path(1,:) + error_x;
% initial_path(1:2, :) = estim_path(1:2, :);

%% -- MANUAL ROTATION -- %%
% NOT USED FOR 06/2023 WORKING VERSION
% v1=0; v2=0; v3=0.2;
% R1=[1 0 0;0 cos(v1) -sin(v1);0 sin(v1) cos(v1)];
% R2=[cos(v2) 0 sin(v2);0 1 0;-sin(v2) 0 cos(v2)];
% R3=[cos(v3) -sin(v3) 0;sin(v3) cos(v3) 0;0 0 1];
% 
% R=R3*R2*R1;
% 
% initial_path=R*initial_path;

% % initial_path_plus(1,:) = estimated_path(1,:)+sigma;
% % initial_path_minus(1,:) = estimated_path(1,:)-sigma;

% model = [true_path, true_path_minus,true_path_plus];

%% -- GENERATE REFERENCE MAP AS ESTIMATE TRAJECTORIES --
model = [estim_path, est_plus_half, est_minus_half, est_plus, est_minus, est_plus2, est_minus2];
% model = [est_plus2];


%% -- NOT CURRENTLY USED --
% NOT USED FOR 06/2023 WORKING VERSION
% data(2,:)=data(2,:)+0.1;
% data(1,:)=data(1,:)+abs(0.2*randn);
% data(2,:)=data(2,:)+abs(0.2*randn);
% data(3,:)=data(3,:)+0.2*randn;

% A plot. Model points and data points in start positions
%%
% NOT USED FOR 06/2023 WORKING VERSION
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


%% -- NOISE -- %%
% NOT USED FOR 06/2023 WORKING VERSION
% fprintf("Initial path: ")
% initial_path(1:3, :)

% for i=1:numel(initial_path(1,:))
%     plus_minus = round(rand);
%     if plus_minus == 0
%         initial_path(3, i) = initial_path(3, i) - rand*5;
%     else
%         initial_path(3, i) = initial_path(3, i) + rand*5;
%     end
% end

%% -- GAUSSIAN NOISE -- 
% NOT USED FOR 06/2023 WORKING VERSION
% sigma_nt = 2;
% N = numel(initial_path(1,:));
% % % noise = sigma_nt*randn(1,N);
% % % for i=1:N
% % %         initial_path(3, i) = initial_path(3, i) + noise(1, i);
% % % end
% % 
% % % sigma_nt = 100;
% % % N = numel(initial_path(1,:));
% noise = sigma_nt*randn(1,N)
% initial_path(3, :) = initial_path(3, :) + noise;
% 
% figure(2)

% fprintf("Initial path with noise: ")
% initial_path(3, :)


%% -- Print Residuals -- %%
% NOT USED FOR 06/2023 WORKING VERSION
% [RotMat,TransVec,dataOut,res]=icp(model,initial_path,[], [], 1,0.00001,true_path);

%% -- CALL ICP FUNCTION -- %%
[RotMat,TransVec,dataOut,res]=icp(model,initial_path,[], [], 1, 0.001, true_path, initial_path, estim_path);


%% 
% NOT USED FOR 06/2023 WORKING VERSION
% disp('Residual is: ')
% res
% disp('Final Matrix X, Y: ')
% dataOut(1:2, 1:end)

% residual = dataOut(1:2,:)-data(1:2, :);
% disp('Final End X, Y: ')
% data(1:2, end)

%% Reference:
%
% Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
% Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2


%% -- PLOT: Model points and data points in transformed positions
figure(3)
% REMAINING CODE IN SECTION NOT USED FOR 06/2023 WORKING VERSION
% subplot(3,1,3)
% f1 = plot3(dataOut(1,:),dataOut(2,:),dataOut(3,:),'b*',true_path(1,:),true_path(2,:),true_path(3,:),'k-o', estim_path(1,:),estim_path(2,:),estim_path(3,:),'r.', initial_path(1,:),initial_path(2,:),initial_path(3,:),'c*')  
% view([0,90])
% 
% % model(1,:),model(2,:),model(3,:),'b*', estim_path(1,:),estim_path(2,:),estim_path(3,:),'g.', true_path(1,:),true_path(2,:),true_path(3,:),'k.', 
% title("Trajectory Fitting Result using ICP")
% xlabel("x(m)");
% ylabel("y(m)");
% legend('ICP solution', 'Truth','Estimate','Initial Position')
% legend('boxoff')
% legend('Location','northeast')

% save_figure("ICP_results",traj_name)

%% -- NOT CURRENTLY USED -- 
% NOT USED FOR 06/2023 WORKING VERSION
% model(1,:),model(2,:),model(3,:),'g.', 
% dataOut(1,:),dataOut(2,:),dataOut(3,:),'b.',
% figure(2)
% plot3(model(1,:),model(2,:),model(3,:),'r.',dataOut(1,:),dataOut(2,:),dataOut(3,:),'g.')
% hold on
% % plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
% title('Transformed data points (green) and model points (red)')

%% 
% NOT USED FOR 06/2023 WORKING VERSION
% true_count = 0;
% estim_count = 0;
% 
% for i=1:numel(initial_path(1,:))
%     plus_minus = round(rand);
%     if plus_minus == 0
%         initial_path(3, i) = initial_path(3, i) - rand*20;
%     else
%         initial_path(3, i) = initial_path(3, i) + rand*10;
%     end
% end


