traj_name = "_3d";

x_range = -1.1:0.05:1.2;
y_range = -3.5:0.05:3.5;

figure (103)
% xlin = linspace(min(xHist), max(xHist), 100);
% ylin = linspace(min(yHist), max(yHist), 100);
[X,Y] = meshgrid(x_range, y_range);
% Z = griddata(x,y,z,X,Y,'natural');
% Z = griddata(x,y,z,X,Y,'cubic');

% Initialize 2D array to store Magnetic Reading [height(Z)] values
ms_fake = zeros(numel(y_range),numel(x_range));

% Compute intensities for the x and y ranges.
k = 0;

for j=1:length(x_range)

    for i=1:length(y_range)
        k = k + 1;
        ms_fake(i,j) = aMap(x_range(j), y_range(i));
        % model(1,k) =  x_range(j);
        % model(2,k) =  y_range(i);
        % model(3,k) =  ms_fake(i,j);   

    end
end

Z = griddata(x_range,y_range,ms_fake,X,Y,'v4');
% mesh(X,Y,Z)
surf(X,Y,Z)



% legend('ICP solution', 'Truth','Estimate','Initial Position')
% legend('boxoff')
% colormap(turbo);

% Use for 2D Contour Map
% contourf(X,Y,Z)

title("3D Magnetic Anomaly Map")
xlabel("x(m)");
ylabel("y(m)");
colorbar;

% For 2D map
% daspect([1 2 1])

% legend('Location','northeast')
axis tight; hold on

save_figure("magnav_refmap", traj_name)

