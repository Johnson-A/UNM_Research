function out = gravity_kernel_function(eval_pts, station_pts, gn)
%GRAVITY_KERNEL_FUNCTION Calculate a linear combination of acquisition
%kernels for all stations over a given set of eval_pts.
%   Since there are more eval_pts than station_pts, kernels are calculated
%   at all eval_pts for every station (in other words, the station index is
%   in the outer loop)

rho_min = 1000;
rho_max = 2000;
C_gn = abs(gn * (rho_max - rho_min));
C_gn_inv = 1 ./ C_gn;
num_pts = size(eval_pts, 2);
num_stations = size(station_pts,2);

out = zeros(num_pts, num_stations);
parfor sn = 1:num_stations
    out(:, sn) = C_gn_inv(sn) * acquisition_kernel(station_pts(:, sn), eval_pts)';
end
out = sum(out, 2) * Constants.G;
end

function out = acquisition_kernel(r_n, r)
    delta_r = bsxfun(@minus, r_n, r);
    dist_3_2 = sum(delta_r .^ 2, 1) .^ (3/2);
    out = 1 ./ (delta_r(3, :) .* dist_3_2);
end

% base = bsxfun(@(station, ~) acquisition_kernel(station_pts(:, station), eval_pts)', ...
%     1:num_stations, (1:num_pts)');
% out = Constants.G * sum(bsxfun(@times, base, C_gn_inv'), 2);
%
% out = Constants.G * arrayfun(@(pt) acquisition_kernel(station_pts, eval_pts(:, pt)) * C_gn_inv, 1:num_pts);
%
% function out = acquisition_kernel_old(C_gn, r_n, r)
%     delta_r = r_n - r;
%     out = Constants.G / C_gn / (delta_r(3) * norm(delta_r)^3);
% end