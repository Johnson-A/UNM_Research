function create_gz_function
    syms x1 x2 y1 y2 z1 z2 y z
    r1(y, z) = gz_vertex(x2, y, z) - gz_vertex(x1, y, z);
    r2(z) = r1(y2, z) - r1(y1, z);
    r3 = r2(z2) - r2(z1);
    matlabFunction(r3 * Constants.G, 'vars', [x1 x2 y1 y2 z1 z2], 'File', 'gz', 'Optimize', true);
end

function gz = gz_vertex(x,y,z)
    r = sqrt(x^2 + y^2 + z^2);
    gz = x * log(y + r) + y * log(x + r) - z * atan(x * y / (z * r));
end