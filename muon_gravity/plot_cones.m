function plot_cones
[X, Y, Elev] = BinaryTerrain.read_file('Tunnel_points_20160715.bin', 23363400, 4600, 5079);

fig = figure(); hold on;
set(gca, 'fontsize', 8)
set(fig, 'Color', 'white');

axis equal tight; lighting gouraud;

% Define sub-section over which to view topography
xind = 1+uint64(size(X,2)*0.0):2:uint64(size(X,2)*1.00);
yind = 1+uint64(size(X,1)*0.2):2:uint64(size(X,1)*0.47);

surf(X(yind, xind), Y(yind, xind), Elev(yind, xind), 'EdgeAlpha', 0.15);
scale = 105;

for mu_station = points_by_regexp(build_table, 'MU[0-9][0-9]')
    [CX, CY, CZ] = cylinder([0, 1], 100);
    
    surf(CX * scale + mu_station(1), CY * scale + mu_station(2), CZ * scale + mu_station(3),...
        'FaceAlpha', 0.2, 'EdgeAlpha', 0.5, 'EdgeColor', 'k', 'FaceColor', 'r');
    
    scatter3(mu_station(1), mu_station(2), mu_station(3), 10, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
end

% Figures generated with dpi > 1200 to avoid graphical artifacts.
dpi = '-r600';

export_fig(fig, 'top_view', '-png', dpi);

view(-100, 6);
export_fig(fig, 'under_view', '-png', dpi);

view(-5, 30);
zticklabels('');
export_fig(fig, 'side_view', '-png', dpi);
end

