#
import csv
import numpy as np
# import magpf 
# import collected trajectories

partial_update = 1

## -- UNCOMMENT THIS SECTION FOR TRAJECTORY FROM FILE -- ##
# TODO Update import file name
traj_name = "trajectory_2.csv"
traj_file = open(traj_name, "r")
trajectory = list(csv.reader(traj_file, delimiter=","))
traj_file.close()

## -- Generate model points to represent Magnetic Map -- ##

## Generate model points to represent Magnetic Map
# Red Points generate from Fake Map (visualize_map.m file)
# Specify the range for plot (size of flight space)
# TODO check, may need different function, linspace maybe?
x_range = [x for x in range(-1.1, 1.2, 0.05)]
y_range = [x for x in range(-3.5, 3.5, 0.05)]

# TODO Update data to be the list of trajectory data
data = trajectory.copy()


# ----- ARTIFICIAL DATA ----- #
sigma = 0.05
true_path = data
initial_path = true_path

# -- GENERATE ESTIMATE PATH AS +SIGMA FROM TRUE PATH --
# WE MAY NEED A MORE SOPHISTICATED APPROACH THAT ACCOUNTS
# FOR CORRELATIONS. MAY BE EXTRACTED FROM COVARIANCE MATRIX 
# OR NORMAL VECTOR TO TRAJECTORY
estim_path = true_path
N = len(estim_path)
n_slopes = N-1
    # Perpendicular Slope
    # Theta
    # Unit Vector cos
    # Unit Vector sin
# TODO check correct conversion of np.zeros
path_data = np.zeros(4,N)


## We must update our true path to have map values.
# Here we update the true_path to have clean mag data
for i in range(N):
    # Overwrite mag measurement
    # TODO need aMap equivalent function for hardware
    true_path[i][2] = aMap(true_path[i][1], true_path[i][2])
    
    
    
## -- CALCULATE SLOPES FROM TRUTH -- 

y_diff = []
x_diff = []

for i in range(n_slopes):
    # Calculate Perpendicular Slope
    
    y_diff.append(estim_path[i+1][2] - estim_path[i][2])
    x_diff.append(estim_path[i+1][1] - estim_path[i][1])
    slope = y_diff / x_diff
    path_data[i][0] = -1 / slope
    # Calculate Theta
    theta = np.atan(path_data[i][0])
    path_data[i][2] = theta
    # Calculate Unit Vector
    path_data[i][2] = np.cos(theta)
    path_data[i][3] = np.sin(theta)
    
path_data[N][3] = path_data[n_slopes][3]
path_data[N][4] = path_data[n_slopes][4]


## CALCULATE ESTIMATED PATH AS OFFSET OF TRUTH
# TODO resume from line 121 of kristy_mag_icp.m