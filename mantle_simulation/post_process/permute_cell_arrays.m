function out = permute_cell_arrays(varargin)
    out = [];
    sizes = cellfun(@(x) 1:numel(x), varargin, 'UniformOutput', false);
    indices = combvec(sizes{:});
    for permutation = indices
        comb = arrayfun(@(ind) varargin{ind}{permutation(ind)}, 1:numel(varargin), 'UniformOutput', false);
        out = [out, comb'];
    end
end