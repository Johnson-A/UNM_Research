function [points, X, Y, Z] = fill_plane(center, x1, x2, extents, n)
%FILL_PLANE generate a set of points in a plane centered at CENTER
%described by the two basis vectors x1, x2. Points are uniformly spaced
%along each basis vector xi with n points between [-EXTENTS(i), +EXTENTS(i)].

    basis_norms = map_columns(@norm, [x1 x2]);
    assert(all(ismembertol([basis_norms - 1, x1' * x2], 0)), ...
           'x1 and x2 must be orthonormal basis vectors');
    
    assert(all(extents > 0), 'Extents must be greater than 0');
    assert(numel(extents) == 2, 'Extents must be specified for both basis vectors');
    
    [XG1, XG2] = meshgrid(linspace(-extents(1), extents(1), n), ...
                          linspace(-extents(2), extents(2), n));
    
    offsets = bsxfun(@times, x1, XG1(:)') + bsxfun(@times, x2, XG2(:)');
    points = bsxfun(@plus, center, offsets);
    reshaped = map_columns(@(col) reshape(col, n, n), points', false);
    [X, Y, Z] = reshaped{:};
end