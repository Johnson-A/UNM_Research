function render_prism(corner, diagonal, xh, yh, zh)
%RENDER_PRISM Render a right rectangular prism in the standard basis.
%   The corner of the prism is given in the standard basis while xh, yh, zh
%   are the relative basis vectors for the oriented prism. The diagonal
%   components are given in the relative basis.

unit_cube = [0 0 0; 0 1 0; 1 1 0; 1 0 0;  % Bottom level, clockwise
             0 0 1; 0 1 1; 1 1 1; 1 0 1]; % Top level, clockwise

scaled_prism = unit_cube .* repmat(diagonal, 1, 8)';
oriented_prism = scaled_prism * [xh yh zh]' + repmat(corner, 1, 8)';

faces = [1 2 3 4;  % 1 |   2
         5 6 7 8;  % 2 |   5
         3 4 8 7;  % 3 | 4 1 3
         1 2 6 5;  % 4 |   6
         2 3 7 6;  % 5 |
         1 4 8 5]; % 6 |


patch('Faces', faces, 'Vertices', oriented_prism, 'FaceColor', 'red', 'EdgeColor', 'w', ...
    'FaceAlpha', 0.1);
end

% render_prism([0;0;0], [1;1;1], [1;1;0] / sqrt(2), [-1;1;0] / sqrt(2), [0;0;-1])