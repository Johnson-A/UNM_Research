function StreamlineTplot_2d
%STREAMLINETPLOT_2D Plot 2D results for mantle simulations.
%   Plot temperature profile of LAB for a given time series. Output the
%   cumulative surface distribution of melt migration along with a
%   comparison of the effects of heat transfer between the solid and fluid
%   phases driven by advection. All physical values used are in SI units.

set(0, 'defaultfigurecolor', 'w');
colormap(parula(1024*16));

output_interval = 5;
start_index = 1;
end_index = 41;
time_range = start_index:output_interval:end_index;

% Root directory containing all data
root_dir = '~/Desktop/mantle_simulation_final_output/'

% All parameter combinations which will be processed
mu_vals = {'5e+21'};
Tbs = {'1300.0'};
k_s = {'0.0', '0.001', '0.01', '0.02'}; % k_s{1} must be 0
k_range = 1:numel(k_s);
assert(isequal(k_s, unique(k_s)));

% Define constants
rho_0 = 3300.0; % SI
alpha = 2.5e-5; % thermal expansion, SI
g     = 9.81;   % SI
kappa_0 = 1.E-6;

% from MultipleRuns.py or codes like it in the folder(s) above, we
% establish the Temp scale
% temp_values = [27.+273, Tb+273, 1300.+273, 1500.+273]
% dTemp = temp_values[3] - temp_values[0]
T_scale = 1305 - 27;
LAB_isotherm = 1573;
h = 1e3; % box dimension in km

hscale = 1e6; % box scale in m
tcont  = 300:100:1600; % in Kelvin
pscale = 1192135725.0; % pressure scale from MultipleRuns.py in Pa

% for streamline calculation, use the following from paraview:
% dimgradP = u*1192135725.0/1e6
% wvel = -(dimgradP-150*9.8*jHat)*1e-15/1e-2
% all SI units

rho_melt = 2800; % kg/m^3

k_over_mu = 1e-13 / 1;

start_x = 1:1:990; % note here units must be in km as displayed in box
start_y = repmat(200, size(start_x));

% nstream = length(start_x) * 0.5;
% seg = 100:50:900;
% startx = repmat(seg, 1, length(seg))';
% nstream = length(startx) * 0.5;
% starty  = 500 * ones(length(startx), 1);
% starty  = reshape(repmat(seg, length(seg), 1), length(seg)^2, 1);
% startz  = 100 * ones(length(startx), 1);
prevDist = [];
numt = 1;

combTracers  = [];
num_bins = 30;
surface_points = [];

    function compare_by_k(do_work, included, fig_base, next_step)
        for ind = 1:numel(included)
            figure(fig_base + ind);
            if exist('next_step', 'var') && next_step, clf; hold on; end
            do_work(included(ind));
        end
    end

    function compare_vs(F1, F2, analysis, c_limits, num_contours, k)
        F_contours = linspace(c_limits(1), c_limits(2), num_contours);
        contourf(X, Y, analysis(F1, F2), F_contours);
        view(2); caxis(c_limits); colorbar; axis equal;
%         title(k);
    end

    function plot_streamlines(sl)
        sl_plot = streamline(sl(1:15:end));

        set(sl_plot, 'color', 'w', 'linewidth', 1.0);
    end

    function plot_surface_points(sp, format)
        sp = sp(end,:);
        bin_dist = linspace(min(X(:)), max(X(:)), num_bins);
        counts = histcounts(sp(:), bin_dist, 'Normalization', 'pdf');
        bin_width = bin_dist(2) - bin_dist(1);
        plot(bin_dist(1:end-1) + bin_width / 2, counts, format);
%         dif = smooth(diff(sp), 15);
%         plot(sp(1:end-1) + (sp(2) - sp(1)) / 2, dif, '-o');
        xlim([bin_dist(1), bin_dist(end)]);
    end

    function sl = calc_streamlines(V)
        VX = V(:,:,1);
        VY = V(:,:,2);

        sl = stream2(X', Y', VX', VY', start_x, start_y, 0.01);
    end

    function setup_single_colorbar(range)
        top = get(subplot(numel(range),1,1), 'Position');
        bot = get(subplot(numel(range),1,numel(range)), 'Position');
        pos_left = 0.5 + (top(4) * 10/4)/2;
        pos_bottom = bot(2);
        width = 0.025;
        height = top(2) + top(4) - pos_bottom;
        colorbar('Position', [pos_left, pos_bottom, width, height]);
    end

    function save_to(fn, step, k, figure_dpi)
        if ~exist('figure_dpi', 'var'), figure_dpi = '-r300'; end

        export_fig(['figures/' fn '_' step 'Myr_k=' k], '-png', figure_dpi);
    end

for vals = permute_cell_arrays(mu_vals, Tbs)
    parsed_vals = num2cell(cellfun(@str2num, vals));
    [mu_scale, ~] = parsed_vals{:};

    base = @(k) [root_dir 'mu=' vals{1} '/Tb=' vals{2} '/k=' k];

    [X,Y,shape] = read_coordinates(base(k_s{1}), h);

    for time_step = time_range
        step_str = num2str(time_step);
        data_set = ['/VisualisationVector/' step_str];

        read_h5 = @(base, fn) h5read([base '/' fn '.h5'], data_set);
        read_into = @(select, dim) @(fn, scale) @(k) ...
            scale * reshape(select(read_h5(base(k), fn)), dim);
        read_scalar = read_into(@(v) v, shape);
        read_vector = read_into(@(v) v', [shape, 3]);

        map_over = @(f, vals) cellfun(f, vals, 'UniformOutput', false);

        M = map_over(read_scalar('mu', mu_scale), k_s);
        T = map_over(read_scalar('T_solid', T_scale), k_s);
        V = map_over(read_vector('v_melt', 1), k_s);

        SL = map_over(@calc_streamlines, V);

        new_streams = map_over(@(sl) cellfun(@(line) line(end, 1), sl), SL);
        surface_points = map_columns(@cell2mat, [surface_points; new_streams], false);

        if time_step == 1, M0 = M; T0 = T; end

        % TODO: Optimize color limits for last step

        for ki = k_range
            'Fractional Temperature Change vs t=0';
            figure(ki + 10); clf; hold on;
            compare_vs(T{ki}, T0{ki}, @(t1, t2) (t1 - t2) ./ t2, [-1, 1] * 0.15, 20, k_s{ki});
            plot_streamlines(SL{ki});
            save_to('fractional_temperature/out', step_str, k_s{ki});
        end

        for ki = k_range(2:end)
            'Temperature difference from k=0 (Degrees C) ';
            figure(ki + 20); clf; hold on;
            compare_vs(T{ki}, T{1}, @(t1, t2) t1 - t2, [-1, 1] * 50, 20, k_s{ki})
            plot_streamlines(SL{ki});
            save_to('delta_temperature/out', step_str, k_s{ki});
        end

        for ki = k_range
            'Fractional Viscosity Change vs t=0';
            figure(ki + 30); clf; hold on;
            compare_vs(M{ki}, M0{ki}, @(t1, t2) log(t1 / t2), [-1, 1] * 1, 20, k_s{ki});
            plot_streamlines(SL{ki});
            save_to('fractional_viscosity/out', step_str, k_s{ki});
        end

        for ki = k_range(2:end)
            'Viscosity difference from k=0 (Degrees C) ';
            figure(ki + 40); clf; hold on;
            compare_vs(M{ki}, M{1}, @(t1, t2) log(t1 / t2), [-1, 1] * 0.8, 20, k_s{ki})
            plot_streamlines(SL{ki});
            save_to('delta_viscosity/out', step_str, k_s{ki});
        end

        %         generate backwards to optimize caxis
        %         if time_step == 1, setup_single_colorbar(c_axis); end

        figure(50); clf; hold on;
        segments = {'-', 'o', '*', 'x'};
        'Cumulative Streamlines Surface Intersections';

        for ki = k_range, plot_surface_points(surface_points{ki}, segments{ki}); end

        legend(k_s{:}); % k{:} or k
        pbaspect([3,1,1]);
        save_to('surface_points/out', step_str, 'all');
    end
end
end

function newTracers = trackStream(newStreams, n, lower, upper)
numStreams = length(newStreams);
xEndPoints = zeros(1,numStreams);

for index = 1:numStreams
    xData = get(newStreams(index), 'XData');
    xEndPoints(index) = xData(end);
end

figure(5); clf; hold on;
subplot(2,1,1); hold on; title('Bin Dist at Current step');
% hist(xEndPoints, n); xlim([lower,upper]);
hist(xEndPoints, 30);

newTracers = histc(xEndPoints, linspace(lower,upper,30));
end

function [X, Y, shape] = read_coordinates(base, scale)
    coords = scale * h5read([base '/T_solid.h5'], '/Mesh/0/coordinates');
    x = coords(1,:);
    y = coords(2,:);

    x_stride = find(x(2:end) < x(1:end-1), 1);

    shape = [x_stride, length(x) / x_stride];

    X = reshape(x, shape);
    Y = reshape(y, shape);
end

% TODO: Alternate implementation, simplified above
%         read_h5 = @(base, fn) h5read([base '/' fn '.h5'], data_set);
%         read_into = @(select, dim) @(fn, scale) @(k) ...
%             scale * reshape(select(read_h5(base(k), fn)), dim);
%         read_scalar = read_into(@(v) v, shape);
%         read_vector = read_into(@(v) v', [shape, 3]);
%
%         map_over = @(f) cellfun(f, k_s, 'UniformOutput', false);
%         to_dict = @(f) containers.Map(k_s, map_over(f));
%
%         M = to_dict(read_scalar('mu'     , mu_scale));
%         T = to_dict(read_scalar('T_solid', T_scale));
%         V = to_dict(read_vector('v_melt' , 1));
%         SL = to_dict(@(k) calc_streamlines(V(k)));
%         new_surface_points = map_over(@(k) cellfun(@(line) line(end,1), SL(k)));
%         surface_points = [surface_points; new_surface_points];
%         if time_step == 1, M0 = M; T0 = T; end
%
%         map_over(@(k) calc_streamlines(V(k), k));
%
%         plot_streamlines = @(k) plot_streamlines(streamlines(k));
%
%         plotT_vs_init = @(k) compare_vs_init(X, Y, T(k), T0(k), k);
%
%         figure_dpi = '-r300';
%
%         figure(1); clf; hold on;
%         title('Temperature Change (Degrees C)')
%         compare_by_k(plotT_vs_init   , k_s);
% %         compare_by_k(plot_streamlines, k_s);
%         export_fig(['figures/fractional_temperature' num2str(time_step)], '-png', figure_dpi);
%
%         figure(2); clf; hold on;
%         title('Temperature Evolution Compared with No Advection')
%         dT_work = @(k) compare_temperature_delta(X, Y, T(k), T('0.0'), k);
%         compare_by_k(dT_work, k_s(2:end));
% %         compare_by_k(plot_streamlines, k_compared);
%         export_fig(['figures/deltaT' num2str(time_step)], '-png', figure_dpi);
%
%         figure(3); clf; hold on;
%
%         map_over(@(k) calc_streamlines(V(k))),
%         title('Streamlines');
%
%         export_fig(['figures/streamlines' num2str(time_step)], '-png', figure_dpi);
