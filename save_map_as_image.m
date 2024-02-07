function save_map_as_image()

%% Saving the figure based on value of beta.
%This block decides if the name of the file is full or partial.
path_to_save = '/home/watersk/mag_ws/map_files/';
name = "mag_map";
traj = "18"
% path_to_save = '';
fig=gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 22 18];
% fig.PaperPosition = [0 0 22 18];
plot_name = sprintf('%s%s_%s.eps',path_to_save,name,traj);
% print('-dpng','-r300',plot_name);      %   *// 300 dpiend
print('-depsc','-r300',plot_name);      %   *// 300 dpiend


end
