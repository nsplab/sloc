costs = dlmread('~/head_corrected_rhs/head.cost_at_grid_pts.10_20.noise_case8.29',' ');
costs_mtx = reshape(costs(:,7), 11,11,11);
slice([-2:0.4:2],[-2:0.4:2],[-2:0.4:2],costs_mtx, [0],[0],[0]);
grid on;
colormap (flipud(jet(64)));
colorbar('vertical');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
daspect([1 1 1]);

figure;
[xi,yi,zi] = meshgrid(-2:.02:2, -2:.02:2, -2:.02:2);
inter_costs_mtx = interp3([-2:0.4:2],[-2:0.4:2],[-2:0.4:2],costs_mtx,xi,yi,zi,'cubic');
slice([-2:0.02:2],[-2:0.02:2],[-2:0.02:2],inter_costs_mtx, [0],[0],[0]);
grid on;
colormap (flipud(jet(64)));
colorbar('vertical');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
daspect([1 1 1]);
shading flat;