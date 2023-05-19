% global partial_update;

function R_yaw = constrainRotation(R)
global partial_update;
% clear all
% close all
% 
% v1 = 0.1;
% v2 = 0.2;
% v3 = 0.3;
% 
% R1=[1 0 0;
%     0 cos(v1) -sin(v1);
%     0 sin(v1) cos(v1)];
% 
% R2=[cos(v2) 0 sin(v2);
%     0 1 0;
%     -sin(v2) 0 cos(v2)];
% 
% R3=[cos(v3) -sin(v3) 0;
%     sin(v3) cos(v3) 0;
%     0 0 1];
% 
% R=R3*R2*R1;

v3 = 0.5*atan(R(2,1)/R(1,1));
% v3 = partial_update*atan2(R(2,1),R(1,1));
% v3 = 0;

R_yaw = [cos(v3) -sin(v3) 0;
    sin(v3) cos(v3) 0;
    0 0 1];



end