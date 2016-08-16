function forward_calculation(n, z0, enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% test_rrpa
% return

% Center prisms
feet_to_meters = 0.3048;
lc = [495740.471; 540894.252; 2116.173];
lc_diag = [12.0; 233.0; 14.0] * feet_to_meters;

large_first_room = lc + [6 - (30 + 8 / 12) / 2; 233; 0] * feet_to_meters;
large_first_room_diag = [30 + 8 / 12, 27 * (1 + 5.5 / 2.25), 16 + 3 / 12] * feet_to_meters;

all_pts = [
    495745.021      540866.046      2115.782;
    495648.807      540807.971      2113.417;
    495748.623      540810.148      2110.525;
    495737.066      540962.890      2116.825;
    495737.109      540952.377      2116.766;
    495726.981      540950.423      2116.765;
    495871.677      541001.772      2215.903;
    495752.092      540864.786      2116.136;
    495606.499      541105.875      2222.672;
    495568.754      541130.415      2223.868;
    495723.007      541022.010      2217.205;
    495723.359      541011.180      2217.036;
    495710.942      541002.394      2214.071;
    495698.732      540998.396      2212.750;
    495690.268      540992.908      2210.663;
    495719.188      540986.871      2207.866;
    495745.642      541071.570      2218.073;
    495738.790      541048.214      2217.552;
    495686.967      541027.777      2217.184;
    495668.580      541083.167      2219.198;
    495645.717      541045.463      2218.234;
    495638.030      541080.819      2219.607;
    495714.539      541069.969      2218.020;
    495729.104      541070.309      2218.092;
    495728.957      541048.231      2217.609;
    495707.473      541135.119      2221.881;
    495845.997      540987.356      2213.440;
    495875.871      540983.127      2213.135;
    495914.670      540973.949      2212.207;
    495951.740      540965.847      2210.869;
    495988.720      540964.036      2209.596;
    495818.250      541043.325      2216.592;
    495751.122      541026.879      2217.016;
    495800.652      541007.167      2214.959;
    495800.649      541007.160      2214.959;
    495677.067      540900.064      2124.393;
    495700.945      540910.649      2135.731;
    495725.015      540920.219      2142.926;
    495743.127      540918.531      2137.259;
    495936.379      541054.789      2215.629;
    495902.716      541109.253      2217.184;
    495872.619      541056.858      2216.713;
    495801.714      541079.070      2217.570;
    495773.663      541119.912      2220.908;
    495738.803      540950.461      2116.762;
    495737.490      540968.022      2116.805;
    495741.057      540921.927      2116.462;
    495743.152      540895.400      2116.205;
    495742.853      540886.973      2116.098;
    495742.279      540894.437      2116.194;
    495741.675      540901.935      2116.254;
    495741.062      540909.402      2116.329;
    495740.503      540916.868      2116.395;
    495740.106      540921.843      2116.465;
    495740.015      540923.108      2116.470;
    495739.927      540924.350      2116.483;
    495739.823      540925.604      2116.483;
    495739.716      540926.831      2116.500;
    495739.532      540929.328      2116.515;
    495739.318      540931.827      2116.541;
    495739.213      540933.075      2116.554;
    495739.105      540934.310      2116.566;
    495738.916      540936.810      2116.606;
    495738.738      540939.284      2116.641;
    495738.553      540941.794      2116.666;
    495738.359      540944.288      2116.687;
    495737.755      540951.733      2116.767;
    495737.207      540959.218      2116.814;
    495736.854      540966.716      2116.802;
    495737.045      540961.244      2116.819;
    495738.871      540961.130      2116.825;
    495735.246      540960.846      2116.837;
    495742.649      540890.400      2116.102;
    495740.471      540894.252      2116.173;
    495744.124      540894.537      2116.192]';

tunnel_pts = [
    495742.853	540886.973	2116.098;
    495742.279	540894.437	2116.194;
    495741.675	540901.935	2116.254;
    495741.062	540909.402	2116.329;
    495740.503	540916.868	2116.395;
    495740.106	540921.843	2116.465;
    495740.015	540923.108	2116.470;
    495739.927	540924.350	2116.483;
    495739.823	540925.604	2116.483;
    495739.716	540926.831	2116.500;
    495739.532	540929.328	2116.515;
    495739.318	540931.827	2116.541;
    495739.213	540933.075	2116.554;
    495739.105	540934.310	2116.566;
    495738.916	540936.810	2116.606;
    495738.738	540939.284	2116.641;
    495738.553	540941.794	2116.666;
    495738.359	540944.288	2116.687;
    495737.755	540951.733	2116.767;
    495737.207	540959.218	2116.814;
    495736.854	540966.716	2116.802]';

start_y = tunnel_pts(2,1);

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

eval_height = z0 + (YI - start_y) * 0.00;

eval_pts = [XI(:)'; YI(:)'; eval_height(:)'];

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

gz_vals = reshape(gz_vals, [n, n]);
R_earth = 6.371E6;
%TODO: Check to make sure sign is correct
slope_correction = 6.67E-11 * 5.972e24 * (1/R_earth^2 - 1./(R_earth + eval_height).^2);
gz_vals = (gz_vals - slope_correction) * 1E5;

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
scatter3(all_pts(1,:), all_pts(2,:), all_pts(3,:) + 2, 10, 'o', 'MarkerEdgeColor','k',...
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
scatter3(tunnel_pts(1,:), tunnel_pts(2,:), tunnel_pts(3,:));

for i = size(voxel_corners, 2)-1:size(voxel_corners, 2),
    render_prism(voxel_corners(:,i), voxel_diag(:,i));
end
 axis equal;

figure(2);  hold on;
title('Calculated gz (mGal)');
surf(XI, YI, gz_vals, 'EdgeColor', 'none');
contour3(XI, YI, gz_vals, 20, 'k');

scatter3(tunnel_pts(1,:), tunnel_pts(2,:), tunnel_pts(3,:), 'ro');

figure(3); hold on;
title('Calculated vs Measured in tunnel (mGal)');
gravity_vals = interp2(XI, YI, gz_vals, tunnel_pts(1,:), tunnel_pts(2,:));

for i = 1:size(tunnel_pts, 2),
    dist_along_tunnel(i) = norm(tunnel_pts(:,i) - tunnel_pts(:,1));
end

plot(dist_along_tunnel, gravity_vals - max(gravity_vals));

tunnel_gz_vals = -[0, 0.3617577612, 0.6441224628, 0.8932126771, 1.1194476632, 1.2599909842, ...
                 1.3050298357, 1.3357904504, 1.3694712209, 1.401333413, 1.4697930874, ...
                 1.533019724, 1.556593782, 1.5893033535, 1.6552009489, 1.7376961059, ...
                 1.797618551, 1.8611792078, 2.0187802839, 2.1911855435, 2.256744684];

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