%% Particle filter
%Authors: Humberto Ramos
%Created: July 23,2021
function [xtrue,xk_updated,state_errors, Pk_updated,t, trajectory_collect] = pf_robot(map,USE_JOYSTICK);
disp('Particle filter initialized....')
%Setup xbox controller
if strcmp(USE_JOYSTICK,'true')
    disp('Joystic Enabled')
    disp('Left stick: move forward')
    disp('Right stick: rotate');
    disp("Press A to read the magnetometer");
%Create joystick object
joy = HebiJoystick(1);
else 
    buttons(1) = 0;
end

%% Initial conditions for estimator.
x0 = [0;-3.5;pi/2];
% delta-time
DT = 1.0; %Not used for now.
%Covariance definitions.
%In the map 6.9085e+03 is amplitude of the smallest peak.
%The measurement noise std is set to be a percentage of it.
smallest_peak_amplitude = 6.9085e+03;
%Take the noise to be = (factor*smalles peak).
factor = 5/100;
Rk = factor* smallest_peak_amplitude;
%Establish process noise values. These will need tunning.
%std_noises(1,1) is meters of error per meter.
%std_noises(2,2) is meters of error per meter.
%std_noises(3,3) is radians of error per revolution.
std_noises = diag([0.01 0.01 deg2rad(1)]);
Pk = blkdiag(0.05^2, 0.05^2, deg2rad(5)^2 );

%% Dynamics
%We are using delta position approach to avoid time variables.
delta_p = 0.2; %meters
delta_psi = 0; %No rotation for now


% Function handler to discretize dynamics
discrete_dyn = @(x,delta_p,delta_psi) x +[delta_p*cos(x(3));
                                delta_p*sin(x(3));
                                delta_psi] ;

% filter measurement
meas = @(x,v) [x(1);x(2)] + v;

%number of states
n = size(x0,1);
%number of measurements.
nm = size(Rk,1);

%Simulation time.
% Tf = 1000.0;
Tf = 100.0;
N = floor(Tf/DT)+1;
t = linspace(0,Tf,N);

% Preallocate variables that log variables.
xk_updated = zeros(N,n);

% log the propagated values 
xk_propagated = zeros(N,n);

% log the measurements
y_measurement = zeros(N,nm);

% log the true state
xtrue = zeros(N,n);

% posterior covariance
Pk_updated = zeros(N,n,n);

% initialize logs
xtrue(1,:) = x0';
xt = x0;
xkp_m = x0;
xk_updated(1,:) = xt;
Pk_updated(1,:,:) = Pk;
%%Generate particles
nP  = 1000;%number of particles.
xhat = x0';%No noise for now for initial state.
xenk = mvnrnd(xhat, Pk, nP)';% Initialize particles
wk = (1/nP)*ones(1,nP);%Initialize normalized weights
xenk_m = zeros(1,n);%EACH COLUMN IS THE HISTORY OF EACH STATE.

figure (1)
%Initialize legends
hold on
propagated_plot_legend = plot(NaN,NaN,'.w','MarkerSize',5, 'DisplayName','Particles');
true_marker_legend = plot(NaN,NaN,'+r','MarkerSize',10,'LineWidth',5,'DisplayName','True position');
mean_marker_legend = plot(NaN,NaN,'+g','MarkerSize',10,'LineWidth',5,'DisplayName','Estimated position');
measure_marker_legend = plot(NaN,NaN,'+g','MarkerSize',10,'LineWidth',5,'DisplayName','Measurement');
legend([propagated_plot_legend,true_marker_legend,mean_marker_legend])
hold off
grid on

% KRISTY ADDED FOR ICP
icp_measure_count = 0;
icp_measurement_list = [];

trajectory_collect = [];

for k = 2:N
    
    if strcmp(USE_JOYSTICK,'false')
        delta_p = 0.2; %
        delta_psi = 0;
    else
        [axes, buttons, povs] = read(joy);
        delta_p = -0.05*axes(2); %Left stick
        delta_psi = -0.08 * axes(4); %Right stick
        %Measurements are taken with button A (Xbox one controller)
    end
    %Propagate for ONE delta
    [xenk] = propagate_pf(xenk,delta_p,delta_psi,std_noises,[]);
     % simulate the system to take measurement then update.
    xt = discrete_dyn(xt,delta_p,delta_psi);
    xtrue(k,:) = xt';
    %It is recommended to compute the estimates before updating and
    %roughening.
    [xenk_m,Pk] = estimate_mean_cov(xenk,wk);
    %Update plots
    figure (1)
    hold on
    propagated_plot = plot(xenk(1,:),xenk(2,:),'.w','MarkerSize',5, 'DisplayName','Particles','HandleVisibility','off');
    mean_plot = plot(xenk_m(1,:),xenk_m(2,:),'+g','MarkerSize',10,'LineWidth',5,'DisplayName','True','HandleVisibility','off');
    true_marker = plot(xt(1),xt(2),'+r','MarkerSize',10,'LineWidth',5,'DisplayName','True','HandleVisibility','off');
    measure_marker = plot(NaN,NaN,'+g','MarkerSize',10,'LineWidth',5,'DisplayName','Measurement','HandleVisibility','off');
    xlim([-1.2 1.2])
    ylim([-3.5 3.5])
    pause(0.2)
    if strcmp(USE_JOYSTICK,'true')
        [axes, buttons, povs] = read(joy);
    end

    if ((mod(k,3)&& strcmp(USE_JOYSTICK,'false')) || buttons(1)==1)     %Get measurement when pressing A
         disp('Measurement received')
         yk = aMap(xt(1),xt(2)) +  mvnrnd(0,Rk,1)'; %Create noisy measurement.
         y_measurement(k,:) = yk';%Log measurement 

         % Modify this value in order to change noise
         data_pose_and_measurement = [xt(1); xt(2); yk]; %aMap(xt(1),xt(2))  OR yk
         trajectory_collect = [trajectory_collect, data_pose_and_measurement ];
         [wk] = update_pf(xenk,yk,Rk,[],wk,[]);

         % KRISTY ADDED FOR ICP
         % [xt(1); xt(2); yk]
         icp_measurement_list = [icp_measurement_list, data_pose_and_measurement]
         icp_measure_count = icp_measure_count + 1;
         if (icp_measure_count >= 4)
             [x, y] = kristy_mag_icp(icp_measurement_list)
             icp_measure_count = 0;
             icp_measurement_list = [];
         end

         figure (1)
         xlim([-1.2 1.2])
         ylim([-3.5 3.5])
         pause(0.1)
    end
    
    if k~=N %do not erase last plot
        figure (1)
        delete(propagated_plot);
        delete(true_marker);
        delete(measure_marker);
        delete(mean_plot);
    end
  
      %Systematic resampling if too few effective particles.
      %Includes roughening.
      effective_N = 1/sum(wk.^2);
      if effective_N < nP/2
         [xenk,wk]=systematic_resampling(wk,xenk);
      end

      
%     Save state estimate and covariance.       
        xk_updated(k,:) = xenk_m';
        Pk_updated(k,:,:) = Pk;
end
%Save history of state errors
state_errors = xk_updated - xtrue;
end
