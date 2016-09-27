function tunnel_rooms = tunnel_spec(lc, tunnel_angle, tunnel_slope)
%TUNNEL_SPEC full physical specification of tunnel given in NAD83 state
%plane coordinates. These physical dimensions should fully match the tunnel
%slope and orientation.

M = rotate_z(tunnel_angle) * rotate_x(atan(tunnel_slope));
[xh, yh, zh] = deal(M(:,1), M(:,2), M(:,3));
flat_yh = normc(yh .* [1;1;0]);

long_corridor = oriented_prism(lc, [3.66; 70.85; 4.27], xh, yh, zh);

side_corridor_diag = [11.90; 2.75; 4.27];
side_corridor_offset = [-side_corridor_diag(1); 59.60; 0];
side_corridor = oriented_prism(lc + M * side_corridor_offset, side_corridor_diag, xh, yh, zh);

side_room_diag = [6.10; 8.10; 4.95];
side_room_offset = side_corridor_offset - ...
    [side_room_diag(1); (side_room_diag(2) - side_corridor_diag(2)) / 2; 0];
side_room = oriented_prism(lc + M * side_room_offset, side_room_diag, xh, yh, zh);

final_room_diag = [9.14; 8.23 + 19.50; 4.95];
final_room_offset = ...
    [-(final_room_diag(1) - long_corridor.diagonal(1)) / 2; long_corridor.diagonal(2); 0];
final_room = oriented_prism(lc + M * final_room_offset, final_room_diag, xh, flat_yh, [0;0;1]);

tunnel_rooms = [long_corridor, side_corridor, side_room, final_room];
end

function R_x = rotate_x(theta)
R_x = [1, 0         ,  0         ;
       0, cos(theta), -sin(theta);
       0, sin(theta),  cos(theta)];
end

function R_z = rotate_z(theta)
R_z = [cos(theta), -sin(theta), 0;
       sin(theta),  cos(theta), 0;
       0         ,  0         , 1];
end