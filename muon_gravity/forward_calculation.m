function forward_calculation(n, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
close all;
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% [X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88.bin', 2422420, 1540, 1573);
[X, Y, Elev] = BinaryTerrain.read_file('Tunnel_points_20160715.bin', 23363400, 4600, 5079);

% Generate submesh on which to downsample
[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));

% Interpolate linearly since we're downsampling
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

min_z = min(ElevI(:));

rho = repmat(Constants.rock_density, n*n, 1);

% eval_pts = [Constants.base_station, Constants.tunnel_pts];

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';

voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];

voxel_diag = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - min_z];

tic;
interaction_matrix = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag);
toc;

lc = point_table{'W wall tunnel', Constants.xyz_index}';

tunnel_rooms = tunnel_spec(lc, Constants.tunnel_angle_offset_from_north, Constants.tunnel_slope);

ind = 1;
for pt = eval_pts
    for p_id = 1:4
        tunnel_effect(ind, p_id) = tunnel_rooms(p_id).eval_gz_at(pt);
    end
    ind = ind + 1;
end

rho_oriented = repmat(-Constants.rock_density, numel(tunnel_rooms), 1);

gz_vals = interaction_matrix * rho + tunnel_effect * rho_oriented;
% inverse = interaction_matrix \ gz_vals;
% diff = sum(abs(inverse - rho)./rho) / numel(rho)

offset_gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;

if ~enable_plotting
    return
end

% Increase the fineness of the default color map
set(0, 'DefaultFigureColormap', parula(1024*16));

elevations = eval_pts(3, :);
northing = eval_pts(2, :);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);

below_cutoff_height = elevations < 2150;

    function do_plot(fig_num, fig_name, mask, N, gz_calc, gz_meas, error)
        figure(fig_num); hold on;
        scatter(N(mask), gz_calc(mask))
        errorbar(N(mask), gz_meas(mask), error(mask), 'o');
        
        title([fig_name ' n = ' num2str(n) ', density = ' num2str(Constants.rock_density)]);
        legend('Calculated', 'Observed');
        xlabel('Northing (m)'); ylabel('gz (mgal)');
        % saveas(gcf, ['figures/' fig_name ' stations_' num2str(n) '_' num2str(int64(Constants.rock_density))], 'png');
    end

do_plot(10, 'Lower',  below_cutoff_height, northing, offset_gz_vals, gz_avg_at_stations, gz_error_at_stations);
do_plot(11, 'Upper', ~below_cutoff_height, northing, offset_gz_vals, gz_avg_at_stations, gz_error_at_stations);

%%
figure(2); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');

surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:) + 2, 10, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');

axis equal tight
lighting gouraud

%%
figure(3); hold on; axis equal;
tunnel_pts = points_by_regexp(point_table, 'TS[0-9][0-9]');

scatter3(tunnel_pts(1,:), tunnel_pts(2,:), tunnel_pts(3,:));

arrayfun(@render, tunnel_rooms);

[p, X, Y, Z] = fill_plane(lc + [0;0;50], [0;1;0], [0;0;1], [300,300], 500);
tic;

inc = 40;
resolution = gravity_kernel_function(p, eval_pts(:,inc), gz_vals(inc) / Constants.rock_density);
toc;
surf(X, Y, Z, log10(reshape(abs(resolution), size(X))), 'EdgeAlpha', 0.1);
scatter3(eval_pts(1,inc), eval_pts(2,inc), eval_pts(3,inc), 10, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
colorbar;
end

function test_rrpa
    x = -150:1:150;
    y = -150:1:150;
    [X, Y] = meshgrid(x,y);

    corner = [-100; -100; -200];
    diagonal = [200; 200; 100];

    m = create_interaction_matrix([X(:)'; Y(:)'; 0 * X(:)'], corner, diagonal);

    calc_gz = reshape(m * 2000, [length(y), length(x)]) * 1E5;

    surf(X, Y, calc_gz, 'EdgeColor', 'none'); hold on;
    contour3(X, Y, calc_gz, 'k');
    axis equal
end