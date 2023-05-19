function [xenk,wk]=systematic_resampling(wk,xenk)
nP = numel(wk);%number of particles.
v = rand();%random number to generate the random positions on the cells.
xenk_new = zeros(size(xenk));
    %Generate the points that will be along the grid of cummulative sum diagram.
    u = zeros(1,nP);
    for j = 1:nP
        %Distributes nP-1 points along the cumulative sum.
        %Note that we add the random offset v to each point.
        u(j) = ((j-1)+v)/nP;
      %For example, with no random v the distribution could look like the
    %following
 %u(1)= 0  u(2)= 1/nP | u(2)= 2/nP | u(3)= 3/nP  u(4)= 4/nP | ...  u(nP-1)= nP-1/nP
   %Now let's see which weights fall into each division
    end
    
    z = cumsum(wk); %Compute cummulative sum of the weights.
    %The cumulative sum will generate numbers that can be used for data
    %segmentation
    
   i =1;
   j = 1;
while j <= nP
    %We keep assigning particles to the ith division, until there is one
    %weight that is outside the division.
   if u(j)< z(i) % z(i) is the ith division in the cumulative sum.
       %This line copies the particle into a new array. It will keep
       %copying the same particle to the new array until no more weights
       %appear in the current ith division.
           xenk_new(:,j) = xenk(:,i);
           j = j + 1;
   else
       %We enter here when we need to move to the next division.
       i = i + 1;
   end
end
    %At this point, xenk_new contains the new cloud of particles.
    %These new set of particles contain copies of the initial particles.
    %The number of copies of each particle, depends on how likely each
    %particle was to represent the state of the system.
    xenk = xenk_new;
   
%      %reset the weights
    wk = (1/nP)*ones(1,nP);
%     
    %Roughening. See Bootstrap Filter from Crassidis and Junkins.
    G = 0.2;
    for i = 1:3
        E(i) = max(xenk(:,i))-min(xenk(:,i));
    end
    cov = (G*E*nP^(-1/3)).^2;
    P_sigmas = diag(cov);
     
    
    for i = 1 :nP
        xenk(:,i) = xenk(:,i) + mvnrnd(zeros(1,3),P_sigmas,1)';
    end
    
    
end
