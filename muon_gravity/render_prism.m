function [ output_args ] = render_prism(corner, diagonal)
%RENDER_PRISM Summary of this function goes here
%   Detailed explanation goes here

unit_cube = [0 0 0;0 1 0;1 1 0;1 0 0;0 0 1;0 1 1;1 1 1;1 0 1];
scaled_prism = unit_cube .* repmat(diagonal, 1, 8)' + repmat(corner, 1, 8)';
faces = [1 2 3 4;5 6 7 8;3 4 8 7;1 2 6 5;2 3 7 6;1 4 8 5];


patch('Faces', faces, 'Vertices', scaled_prism, 'FaceColor', 'red', 'EdgeColor', 'w', ...
    'FaceAlpha', 0.1);
end

