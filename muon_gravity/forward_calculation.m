function forward_calculation(n, z0, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

start_x = 495741.5;
start_y = 540886;

file_name = 'TA41_Tunnel_LIDAR_NAVD88';

[X, Y, Elev] = read_binary_file(file_name);

% Generate submesh on which to downsample
[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));

% Interpolate linearly since we're downsampling
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

min_z = min(ElevI(:));

rho = 2600 * ones(n*n, 1); %  + 2600 * (rand(n*n, 1) - 0.5)

eval_height = z0 + (YI - start_y) * 0.0;

eval_pts = [XI(:)'; YI(:)'; eval_height(:)'];

along_row = ones(1, n*n);
voxel_corners = [XI(:)'; YI(:)'; min_z * along_row];
voxel_diag = [dx * along_row; dy * along_row; ElevI(:)' - min_z];

interaction_matrix = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag);

gz_vals = interaction_matrix * rho;
inverse = interaction_matrix \ gz_vals;
diff = sum(abs(inverse - rho)./rho) / numel(rho)

gz_vals = reshape(gz_vals, [n, n]);
R_earth = 6.371E6;
%TODO: Check to make sure sign is correct
slope_correction = 6.67E-11 * 5.972e24 * (1/R_earth^2 - 1./(R_earth + eval_height).^2);
gz_vals = (gz_vals - slope_correction) * 1E5;

if ~enable_plotting
    return
end

colormap(parula(1024));

figure(1); hold on;
title('Interpolated Elevation Data');
surf(XI, YI, ElevI, 'EdgeAlpha', 0.2);

figure(2);  hold on;
title('Calculated gz (mGal)');
surf(XI, YI, gz_vals, 'EdgeColor', 'none');
contour3(XI, YI, gz_vals, 20, 'k');

tunnel_x = linspace(start_x, start_x, 100);
tunnel_y = linspace(start_y, start_y + 83, 100);

figure(3); hold on;
title('Calculated vs Measured in tunnel (mGal)');
gravity_vals = interp2(XI, YI, gz_vals, tunnel_x, tunnel_y);
plot(tunnel_y - start_y, gravity_vals - max(gravity_vals));

tunnel_x_pts = [0, 7.5, 15, 22.5, 30, 35, 36.25, 37.5, 38.75, 40, 42.5, 45, 46.25, 47.5, 50, 52.5, 55, 57.5, 65, 72.5, 80];
tunnel_y_pts = -[0, 0.3617577612, 0.6441224628, 0.8932126771, 1.1194476632, 1.2599909842, ...
                 1.3050298357, 1.3357904504, 1.3694712209, 1.401333413, 1.4697930874, ...
                 1.533019724, 1.556593782, 1.5893033535, 1.6552009489, 1.7376961059, ...
                 1.797618551, 1.8611792078, 2.0187802839, 2.1911855435, 2.256744684];

scatter(tunnel_x_pts, tunnel_y_pts);
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
    topo = importdata([file_name '.txt'], ',', 1);
    fileID = fopen([file_name '.bin'],'w');
    fwrite(fileID, topo.data, 'single');
    fclose(fileID);
end

function [X,Y,Elev] = read_binary_file(file_name)
    total_points = 2422420;
    num_points_x = 1540;
    num_points_y = 1573;
    assert(num_points_x * num_points_y == total_points);

    fileID = fopen([file_name '.bin']);
    topo = fread(fileID, [total_points, 3], 'single') * 0.3048;
    fclose(fileID);

    X    = reshape(topo(:,2), [num_points_x, num_points_y])';
    Y    = reshape(topo(:,3), [num_points_x, num_points_y])';
    Elev = reshape(topo(:,1), [num_points_x, num_points_y])';
end