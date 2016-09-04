function forward_calculation(n, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
close all;
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% [X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88.bin', 2422420, 1540, 1573);
[X, Y, Elev] = read_binary_file('Tunnel_points_20160715.bin', 23363400, 4600, 5079);

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
xyz_index = {'Easting', 'Northing', 'Elevation'};

eval_pts = point_table{measured_points, xyz_index}';

voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];

voxel_diag = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - min_z];

tic;
interaction_matrix = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag);
toc;

lc = point_table{'W wall tunnel', xyz_index}';

tunnel_rooms = tunnel_spec(lc, Constants.tunnel_angle_offset_from_north, Constants.tunnel_slope);

ind = 1;
for pt = eval_pts
    for p_id = 1:4
        tunnel_effect(ind, p_id) = tunnel_rooms(p_id).eval_gz_at(pt);
    end
    ind = ind + 1;
end

rho_oriented = repmat(-Constants.rock_density, numel(tunnel_rooms), 1);
% rho_oriented = [-Constants.rock_density; -Constants.rock_density + 500];

gz_vals = interaction_matrix * rho + tunnel_effect * rho_oriented;
% inverse = interaction_matrix \ gz_vals;
% diff = sum(abs(inverse - rho)./rho) / numel(rho)

gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;

if ~enable_plotting
    return
end

% Increase the fineness of the default color map
set(0, 'DefaultFigureColormap', parula(1024*16));

elevations = eval_pts(3, :);
northing = eval_pts(2, :);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

% TODO fix plotting code
gz_avg_at_stations = cellfun(@(c) mean(c), measured_values);
gz_error_at_stations = cellfun(@(c) norm(c), measure_errors);

below_cutoff_height = elevations < 2150;

figure(10); hold on;
scatter(northing(below_cutoff_height), gz_vals(below_cutoff_height))
errorbar(northing(below_cutoff_height), gz_avg_at_stations(below_cutoff_height), ...
    gz_error_at_stations(below_cutoff_height), 'o');

title(['n = ' num2str(n) ', density = ' num2str(Constants.rock_density)]);
legend('Calculated', 'Observed');
xlabel('Northing (m)'); ylabel('gz (mgal)');
% saveas(gcf, ['figures/Lower stations_' num2str(n) '_' num2str(int64(Constants.rock_density))], 'png');

figure(11); hold on;
scatter(northing(~below_cutoff_height), gz_vals(~below_cutoff_height))
errorbar(northing(~below_cutoff_height), gz_avg_at_stations(~below_cutoff_height), ...
    gz_error_at_stations(~below_cutoff_height), 'o');

title(['n = ' num2str(n) ', density = ' num2str(Constants.rock_density)]);
legend('Calculated', 'Observed');
xlabel('Northing (m)'); ylabel('gz (mgal)');
% saveas(gcf, ['figures/Upper stations_' num2str(n) '_' num2str(int64(Constants.rock_density))], 'png');

figure(2); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');
% for i = 1:size(voxel_corners, 2),
%     render_prism(voxel_corners(:,i), voxel_diag(:,i), [1;0;0], [0;1;0], [0;0;1]);
% end

surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:) + 2, 10, 'o', ...
    'MarkerEdgeColor','k', 'MarkerFaceColor', 'r');

axis equal tight
lighting gouraud

figure(3); hold on; axis equal;
is_tunnel_pt = ~cellfun(@isempty, regexp(point_table.Properties.RowNames, 'TS[0-9][0-9]'));
tunnel_pts = point_table{is_tunnel_pt, xyz_index}';

scatter3(tunnel_pts(1,:), tunnel_pts(2,:), tunnel_pts(3,:));

for prism = tunnel_rooms
    prism.render
end
% saveas(gcf, 'figures/Tunnel_Spec', 'png');

% figure(2);  hold on;
% title('Calculated gz (mGal)');
% surf(XI, YI, gz_vals, 'EdgeColor', 'none');
% contour3(XI, YI, gz_vals, 20, 'k');
%
% scatter3(Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:), Constants.tunnel_pts(3,:), 'ro');

% figure(3); hold on;
% title('Calculated vs Measured in tunnel (mGal)');
% gravity_vals = interp2(XI, YI, gz_vals, Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:));

% for i = 1:size(Constants.tunnel_pts, 2),
%     dist_along_tunnel(i) = norm(Constants.tunnel_pts(:,i) - Constants.tunnel_pts(:,1));
% end
%
% plot(dist_along_tunnel, gravity_vals - max(gravity_vals));
%
% scatter(dist_along_tunnel, tunnel_gz_vals);
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

function write_binary_file_from_txt(file_name)
    fid = fopen(file_name, 'r');
    data = textscan(fid, '%d %f %f %f', 'HeaderLines', 1, 'Delimiter', ',');
    fclose(fid);

    binary_file = fopen([file_name '.bin'], 'w');
    fwrite(binary_file, [data{2}, data{3}, data{4}], 'single');
    fclose(binary_file);
end

function [X,Y,Elev] = read_binary_file(file_name, total_points, num_points_x, num_points_y)
    assert(num_points_x * num_points_y == total_points);

    fileID = fopen(file_name);
    topo = fread(fileID, [total_points, 3], 'single') * 0.3048;
    fclose(fileID);

    X    = reshape(topo(:,2), [num_points_x, num_points_y])';
    Y    = reshape(topo(:,3), [num_points_x, num_points_y])';
    Elev = reshape(topo(:,1), [num_points_x, num_points_y])';
end