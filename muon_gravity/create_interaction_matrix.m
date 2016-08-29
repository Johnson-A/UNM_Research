% TODO: Paralellize by whichever there are more of?
% Store 3-vector along column since matlab stores in column-major order
function m = create_interaction_matrix(eval_pts, voxels)
    num_pts = size(eval_pts, 2);
    num_voxels = numel(voxels);
    m = zeros(num_pts, num_voxels);

    parfor voxel_id = 1:num_voxels,
        voxel = voxels(voxel_id);
        
        for pt = 1:num_pts,
            m(pt, voxel_id) = voxel.eval_gz_at(eval_pts(:, pt));
        end
    end
end