function gz = oriented_prism_gz(voxel_corner, voxel_diag, xh, yh, zh, eval_pt)
%ORIENTED_PRISM_GZ GZ calculation for an arbitrarily oriented right
%rectangular prizm
%   The voxel_corner is the same regardless

% Calculate the full gravitational acceleration in the euclidean frame of
% the prism


% Find the eval_pt location in the prism's euclidean orientation
offset_pt = eval_pt - voxel_corner;

A = [xh, yh];
xy = pinv(A) * offset_pt; % Moore Penrose pseudoinverse
z = dot(offset_pt - A * xy, zh);

corner_offset = -[xy; z];

bounds = [corner_offset, corner_offset + voxel_diag];

g = full_g_vector(bounds(1, :), bounds(2, :), bounds(3, :));
z_dir = [0; 0; 1];
z_dir_in_prism_basis = [dot(z_dir, xh); dot(z_dir, yh); dot(z_dir, zh)];
gz = dot(g, z_dir_in_prism_basis);
end

function g = full_g_vector(x, y, z)
g = [gz_mixed(y, z, x);
     gz_mixed(z, x, y);
     gz_mixed(x, y, z)];
end

function gz_val = gz_mixed(xb, yb, zb)
    gz_val = gz(xb(1), xb(2), yb(1), yb(2), zb(1), zb(2));
end

