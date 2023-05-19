%% plot_ellipse.m
%This function plots the error ellipse that corresponds to the error covariance
%matrix Px centered at the mean mu.
function plot_ellipse(Px,mu,str_filter_name,color)
% Calculate the eigenvectors and eigenvalues.
%eigenvectors in columns of V.
[V, D ] = eig(Px);

% First max gets maximum for each colum. Then the maximum of the maximum is
% computed.
[i, j] = find(D == max(max(D)));
%Look for the corresponding eigenvalue.
largest_V = V(:, i);

% Get the largest eigenvalue
largest_eigenval = max(max(D));

% Get the smallest eigenvector and eigenvalue
if(i == 1)
    smallest_eigenval = max(D(:,2));
    smallest_eigenvec = V(:,2);
else
    smallest_eigenval = max(D(:,1));
    smallest_eigenvec = V(1,:);
end

%Ellipse orientation.
angle = atan2(largest_V(2), largest_V(1));

%Keep ellipse orientation between 0 and 2pi
if(angle < 0)
    angle = angle + 2*pi;
end

% Get the 95% confidence interval error ellipse
chisquare_val = 2.4477;
theta_grid = linspace(0,2*pi);
phi = angle;
X0=mu(1);
Y0=mu(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Rotation matrix 
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%Rotate coordinates
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
hold on;
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,color)
mark = sprintf(strcat(color,'+'));
plot (X0,Y0,mark,'LineWidth',3)

legend(sprintf('Ellipse error and estimated value %s',str_filter_name))

end