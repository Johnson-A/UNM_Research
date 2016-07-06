% TODO: Paralellize by whichever there are more of?
% Store 3-vector along column since matlab stores in column-major order
function m = create_interaction_matrix(eval_pts, voxel_corner, voxel_diag)
    num_pts = size(eval_pts, 2);
    num_voxels = size(voxel_corner, 2);
    m = zeros(num_pts, num_voxels);

    parfor voxel_id = 1:num_voxels,
        corner = voxel_corner(:, voxel_id);
        diag = voxel_diag(:, voxel_id);
        
        for pt = 1:num_pts,
            c = corner - eval_pts(:, pt);

            m(pt, voxel_id) = gz(c(1), c(1) + diag(1), ...
                                 c(2), c(2) + diag(2), ...
                                 c(3), c(3) + diag(3)) * 6.67E-11;
        end
    end
end