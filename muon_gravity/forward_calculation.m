function forward_calculation(n, z0)
%FORWARD_CALCULATION Terrain based forward model gravity calculation

% test_rrpa
% write_binary_file_from_txt

colormap(parula(1024));
file_name = 'TA41_Tunnel_LIDAR_NAVD88';

[X,Y,Elev] = read_binary_file(file_name);
% surf(X,Y,Elev,'EdgeAlpha', '0.2')

% Resize using the default nearest neighbor algorithm
[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));

% Interpolate linearly
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');

surf(XI, YI, ElevI, 'EdgeAlpha', 0.2);

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

xpts = reshape(XI, [1, n * n]);
ypts = reshape(YI, [1, n * n]);
elev_vec = reshape(ElevI, [1, n * n]);  
min_z = min(ElevI(:));

rho = 2600 * ones(n*n, 1);

interaction_matrix = zeros(n*n);

grid = [xpts; ypts];
num_pts = length(grid);
eval_height = z0 + 0.00001 + (YI - 540886) * 0.01;

eval_height_vec = reshape(eval_height, [1, n*n]);

syms x1 x2 y1 y2 z1 z2 y z
r1(y, z) = c_fun(x2, y, z) - c_fun(x1, y, z);
r2(z) = r1(y2, z) - r1(y1, z);
r3 = 6.67E-11 * (r2(z2) - r2(z1));
gz = matlabFunction(r3, 'vars', [x1 x2 y1 y2 z1 z2]);

parfor voxel_id = 1:num_pts,
    corner = [grid(:, voxel_id); min_z];
    diagonal = [dx; dy; elev_vec(voxel_id) - min_z];
    
    for grid_pt = 1:num_pts,
        p = [grid(:, grid_pt); eval_height(grid_pt)];
        
%         contribution = right_rectangular_prism_acceleration(corner, diagonal, p);
%         contribution = rewrite(corner - p, diagonal);
%         interaction_matrix(grid_pt, voxel_id) = contribution;
        c1 = corner - p;
        
        interaction_matrix(grid_pt, voxel_id) = gz(c1(1), c1(1) + diagonal(1), ...
                                                   c1(2), c1(2) + diagonal(2), ...
                                                   c1(3), c1(3) + diagonal(3));
    end
end

gz = interaction_matrix * rho;
% test = interaction_matrix \ gz

gz = reshape(gz, [n, n]);
figure();
R_earth = 6.371E6;
slope_correction = 6.67E-11 * 5.972e24 * (1/R_earth^2 - 1./(R_earth + eval_height).^2);
gz = (gz - slope_correction) * 1E5;

surf(XI, YI, gz, 'EdgeColor', 'none'); hold on;
contour3(XI, YI, gz, 20, 'k');

tunnel_x = linspace(495741.5, 495741.5, 100);
tunnel_y = linspace(540886, 540886 + 83, 100);

x_width = max(XI(:)) - min(XI(:));
y_width = max(YI(:)) - min(YI(:));
disp(495741.5 - 2 * x_width);
disp(495741.5 + 2 * x_width);
disp(540886 - 2 * y_width);
disp(540886 + 2 * y_width);

figure(); hold on;
gravity_vals = interp2(XI, YI, gz, tunnel_x, tunnel_y);
plot(tunnel_y - 540886, gravity_vals - max(gravity_vals));

tunnel_x_pts = [0, 7.5, 15, 22.5, 30, 35, 36.25, 37.5, 38.75, 40, 42.5, 45, 46.25, 47.5, 50, 52.5, 55, 57.5, 65, 72.5, 80];
tunnel_y_pts = -[0, 0.3617577612, 0.6441224628, 0.8932126771, 1.1194476632, 1.2599909842, ...
                 1.3050298357, 1.3357904504, 1.3694712209, 1.401333413, 1.4697930874, ...
                 1.533019724, 1.556593782, 1.5893033535, 1.6552009489, 1.7376961059, ...
                 1.797618551, 1.8611792078, 2.0187802839, 2.1911855435, 2.256744684];

scatter(tunnel_x_pts, tunnel_y_pts);
end

% Is calculating inside a prism okay?
function gz = right_rectangular_prism_acceleration(corner, diagonal, p)
    G = 6.67E-11;
    gz = 0;

    for i = [1,2],
        for j = [1,2],
            for k = [1,2],
                u = (-1)^i * (-1)^j * (-1)^k;
                vertex = corner + [(i-1); (j-1); (k-1)] .* diagonal;
                delta = vertex - p;
                dxi = delta(1); dyj = delta(2); dzk = delta(3);
                R_ijk = norm(delta, 2);

                gz = gz + G * u * ( -dzk * atan(dxi * dyj / (dzk * R_ijk)) ...
                    + dxi * log(R_ijk + dyj) ...
                    + dyj * log(R_ijk + dxi));
            end
        end
    end
end

function gz = c_fun(x,y,z)
    r = sqrt(x^2 + y^2 + z^2);
    gz = x * log(y + r) + y * log(x + r) - z * atan(x * y / (z * r));
end

function gz = rewrite(corner, diag)
    r1 = @(y,z) c_fun(corner(1) + diag(1), y, z) - c_fun(corner(1), y, z);
    r2 = @(z) r1(corner(2) + diag(2), z) - r1(corner(2), z);
    gz = 6.67E-11 * (r2(corner(3) + diag(3))- r2(corner(3)));
end

function test_rrpa
    x = -150:1:150;
    y = -150:1:150;
    corner = [-100; -100; -200];
    diagonal = [200; 200; 100];
    res = zeros(length(x), length(y));

    for i = 1:length(x),
        for j = 1:length(y),
            res(i,j) = right_rectangular_prism_acceleration(corner, diagonal, [x(i); y(j); 0]);
        end
    end

    surf(2000 * res * 1E5, 'EdgeColor', 'none'); hold on;
    contour3(2000 * res * 1E5, 'k');
    axis equal
end

function write_binary_file_from_txt
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