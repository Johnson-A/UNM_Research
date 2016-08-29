function forward_calculation(n, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% [X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88', 2422420, 1540, 1573);
[X, Y, Elev] = read_binary_file('Tunnel_points_20160715', 23363400, 4600, 5079);

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

eval_pts = [Constants.base_station, Constants.tunnel_pts];

voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];

voxel_diag = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - min_z];
voxels = repmat(standard_prism([0;0;0], [0;0;0]), n*n, 1);

for ind = 1:n*n,
    voxels(ind) = standard_prism(voxel_corners(:, ind), voxel_diag(:, ind));
end

interaction_matrix = create_interaction_matrix(eval_pts, voxels);

tunnel_effect = create_interaction_matrix(eval_pts, Constants.tunnel_rooms);

rho_oriented = repmat(-Constants.rock_density, 4, 1);
% rho_oriented = [-Constants.rock_density; -Constants.rock_density + 500];

gz_vals = interaction_matrix * rho + tunnel_effect * rho_oriented;
inverse = interaction_matrix \ gz_vals;
diff = sum(abs(inverse - rho)./rho) / numel(rho)

gz_vals = gz_vals * 1E5;

if ~enable_plotting
    return
end

colormap(parula(1024*16));

figure(1);
plot(gz_vals - gz_vals(1), 'o-'); hold on;
tunnel_gz_vals = [
    9.13086E-15
    -0.399272713
    -0.768195747
    -1.048850725
    -1.297977524
    -1.523139457
    -1.668033758
    -1.703311442
    -1.734109862
    -1.767828446
    -1.802984822
    -1.868226903
    -1.941130622
    -1.955196262
    -1.987952838
    -2.068448743
    -2.136413966
    -2.204271225
    -2.259974671
    -2.430427306
    -2.59819765
    -2.669678293];

plot(tunnel_gz_vals, 'o-')
legend('Calculated', 'Observed');
title('Calculated vs Observed gz values');
xlabel('pt'); ylabel('gz (mgal)');
axis tight

figure(2); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');
% for i = 1:size(voxel_corners, 2),
%     render_prism(voxel_corners(:,i), voxel_diag(:,i), [1;0;0], [0;1;0], [0;0;1]);
% end

surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
scatter3(Constants.all_pts(1,:), Constants.all_pts(2,:), Constants.all_pts(3,:) + 2, 10, 'o', 'MarkerEdgeColor','k',...
        'MarkerFaceColor', 'r');

axis equal tight
lighting gouraud

figure(3); hold on; axis equal;
scatter3(Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:), Constants.tunnel_pts(3,:));

for prism = Constants.tunnel_rooms
    prism.render
end

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
    fid = fopen([file_name '.txt'], 'r');
    data = textscan(fid, '%d %f %f %f', 'HeaderLines', 1, 'Delimiter', ',');
    fclose(fid);

    binary_file = fopen([file_name '.bin'], 'w');
    fwrite(binary_file, [data{2}, data{3}, data{4}], 'single');
    fclose(binary_file);
end

function [X,Y,Elev] = read_binary_file(file_name, total_points, num_points_x, num_points_y)
    assert(num_points_x * num_points_y == total_points);

    fileID = fopen([file_name '.bin']);
    topo = fread(fileID, [total_points, 3], 'single') * 0.3048;
    fclose(fileID);

    X    = reshape(topo(:,2), [num_points_x, num_points_y])';
    Y    = reshape(topo(:,3), [num_points_x, num_points_y])';
    Elev = reshape(topo(:,1), [num_points_x, num_points_y])';
end