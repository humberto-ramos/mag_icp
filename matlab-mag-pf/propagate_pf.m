function [xenk]=propagate_pf(xenk,delta_p,delta_psi,Q,DT)
nP = size(xenk,2);
        %Propagate ensemble
        for i = 1: nP
            %This block contains time. We may use this code later on.
%             delta_theta =  omega*DT;
%             xenk(3,i) = xenk(3,i) + delta_theta + (delta_theta/(2*pi))*Q(3,3)*randn();
%             delta_distance = V*DT + V*DT*Q(1,1)*randn();%Nominal distance (V*DT)+  
%             xenk(1:2,i)  =  xenk(1:2,i)+ [ delta_distance*cos(xenk(3,i));
%                                       delta_distance*sin(xenk(3,i))] ;
%           
            %For now we will use delta_p and delta_psi only.
            %CONTINUE HERE
            xenk(3,i) = xenk(3,i) + delta_psi + Q(3,3)*randn();
            delta_distance = delta_p + Q(1,1)*randn();%Nominal distance (V*DT)+  
            xenk(1:2,i)  =  xenk(1:2,i)+ [ delta_distance*cos(xenk(3,i));
                                      delta_distance*sin(xenk(3,i))] ;
                                  
                                  
        end
   
