function magnetometer_reading = meas_from_map(xhat,F)
%This function generates an artificial mag measurement based on the known map.
%INPUT: state vector, map
%OUTPUT: Expected magnetometer reading
%Use: magnetometer_reading = meas_from_map(xhat,map)
%Magnetometer reading in X
magnetometer_reading(1,1) = F.x([xhat(1),xhat(2),xhat(3)]);
% magnetometer_reading(1)=  griddatan(,,);
%Magnetometer reading in Y
magnetometer_reading(2,1) = F.y([xhat(1),xhat(2),xhat(3)]);
%Magnetometer reading in Z
magnetometer_reading(3,1) = F.z([xhat(1),xhat(2),xhat(3)]);

end

