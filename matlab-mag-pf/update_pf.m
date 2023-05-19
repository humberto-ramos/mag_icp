function [wk] = update_pf(xenk,yk,Rk,meas,wk,F)
nP = size(xenk,2);%number of particles

    for i=1:nP
%         y_expected = meas_from_map(xenk(:,i),F);
        y_expected = aMap(xenk(1,i),xenk(2,i));

        %note that the argument for the measurement function does not
        %include sensor noise. The likelihood will take care of it.
        like = -0.5*(yk-y_expected)'*inv(Rk)*(yk-y_expected);
        
        wk(i) = wk(i)*exp(like);
        %This factor sqrt(det((2*pi)*Rk)) is not included in the
        %likelihood, but it does not matter since it can be factored
        %and then cancelled out during the normalization.
     end
    %weights normalization.
    sum_of_weights = sum(wk);
    wk = wk/sum_of_weights;
    
  
end