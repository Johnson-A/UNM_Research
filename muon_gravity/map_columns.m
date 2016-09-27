function result = map_columns(fun, A, is_scalar)
    if nargin == 2
        is_scalar = true;
    elseif nargin ~= 3
        error('Wrong number of arguments')
    end
    
    result = arrayfun(@(ind) fun(A(:, ind)), 1:size(A, 2), 'UniformOutput', is_scalar);
end