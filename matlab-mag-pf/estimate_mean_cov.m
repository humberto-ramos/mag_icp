function [xmean,Pk] = estimate_mean_cov(xenk,wk)
nP = size(xenk,2);%number of particles

%Compute mean and covariance after the weight update.
    sumx = 0;
    for i=1:nP
        sumx = sumx +  wk(i)*xenk(:,i);
    end
    xmean = sumx;
    
    sumP = 0;
    for i = 1: nP
        dev = xenk(:,i)-xmean;
        sumP = sumP + wk(i)*(dev*dev');
    end 
    Pk = sumP;