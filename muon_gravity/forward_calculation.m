function forward_calculation(n, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% Center prisms
lc = [495740.471; 540894.252; 2116.173];
lc_diag = [12.0; 233.0; 14.0] * Constants.feet_to_meters;

large_first_room = lc + [6 - (30 + 8 / 12) / 2; 233; 0] * Constants.feet_to_meters;
large_first_room_diag = [30 + 8 / 12, 27 * (1 + 5.5 / 2.25), 16 + 3 / 12] * Constants.feet_to_meters;

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

rock_density = 1600;
rho = rock_density * ones(n*n, 1);
rho(end+1) = -rock_density;
rho(end+1) = -rock_density + 500;

% eval_pts = [XI(:)'; YI(:)'; eval_height(:)'];
eval_pts = Constants.tunnel_pts;

% Center voxel instead of at corner
along_row = ones(1, n*n);
voxel_corners = [XI(:)'; YI(:)'; min_z * along_row];
voxel_corners(:,end+1) = lc;
voxel_corners(:,end+1) = large_first_room;

voxel_diag = [dx * along_row; dy * along_row; ElevI(:)' - min_z];
voxel_diag(:,end+1) = lc_diag;
voxel_diag(:,end+1) = large_first_room_diag;

interaction_matrix = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag);

gz_vals = interaction_matrix * rho;
inverse = interaction_matrix \ gz_vals;
diff = sum(abs(inverse - rho)./rho) / numel(rho)

% gz_vals = reshape(gz_vals, [n, n]) * 1E5;
gz_vals = gz_vals * 1E5;

plot(gz_vals - max(gz_vals), 'o-'); hold on;
tunnel_gz_vals = -[0, 0.3617577612, 0.6441224628, 0.8932126771, 1.1194476632, 1.2599909842, ...
                 1.3050298357, 1.3357904504, 1.3694712209, 1.401333413, 1.4697930874, ...
                 1.533019724, 1.556593782, 1.5893033535, 1.6552009489, 1.7376961059, ...
                 1.797618551, 1.8611792078, 2.0187802839, 2.1911855435, 2.256744684];

plot(tunnel_gz_vals, 'o-')
if ~enable_plotting
    return
end

colormap(parula(1024*16));

fig=figure(1); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');
% surf(XI, YI, ElevI, 'EdgeAlpha', 0.2);
% for i = 1:size(voxel_corners, 2),
%     render_prism(voxel_corners(:,i), voxel_diag(:,i));
% end

% xind = uint64(size(X,2)*0.2):2:uint64(size(X,1)*0.47);
% yind = 1:2:size(X, 2);

surf(XI, YI, ElevI, 'EdgeAlpha', 0.2);
scatter3(Constants.all_pts(1,:), Constants.all_pts(2,:), Constants.all_pts(3,:) + 2, 10, 'o', 'MarkerEdgeColor','k',...
        'MarkerFaceColor', 'r');

axis equal tight
lighting gouraud

% [CX, CY, CZ] = cylinder([0, 1], 100);
% surf(CX * 105 + tunnel_pts(1, end), CY * 105 + tunnel_pts(2, end), CZ * 105 + tunnel_pts(3, end),...
%     'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, 'EdgeColor', 'k', 'FaceColor', 'g');
%
% set(gca, 'fontsize', 8)
% set(fig, 'Color', 'white');
% export_fig(fig, 'top_view', '-png', '-r1400');
%
% view(-5, 30);
% export_fig(fig, 'side_view', '-png', '-r1400');
%
% view(-100, 6);
% export_fig(fig, 'under_view', '-png', '-r1400');
%
% return

figure(4); hold on;
scatter3(Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:), Constants.tunnel_pts(3,:));

for i = size(voxel_corners, 2)-1:size(voxel_corners, 2),
    render_prism(voxel_corners(:,i), voxel_diag(:,i));
end
axis equal;

figure(2);  hold on;
title('Calculated gz (mGal)');
surf(XI, YI, gz_vals, 'EdgeColor', 'none');
contour3(XI, YI, gz_vals, 20, 'k');

scatter3(Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:), Constants.tunnel_pts(3,:), 'ro');

figure(3); hold on;
title('Calculated vs Measured in tunnel (mGal)');
gravity_vals = interp2(XI, YI, gz_vals, Constants.tunnel_pts(1,:), Constants.tunnel_pts(2,:));

for i = 1:size(Constants.tunnel_pts, 2),
    dist_along_tunnel(i) = norm(Constants.tunnel_pts(:,i) - Constants.tunnel_pts(:,1));
end

plot(dist_along_tunnel, gravity_vals - max(gravity_vals));

scatter(dist_along_tunnel, tunnel_gz_vals);
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