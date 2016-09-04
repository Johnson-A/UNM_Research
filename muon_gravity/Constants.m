classdef Constants
    %CONSTANTS Constant definition file.

    properties (Constant)
        feet_to_meters = 0.3048;
        G = 6.67E-11;
        rock_density = 1.4625E3;

        tunnel_angle_offset_from_north = (4 + 31 / 60 + 27 / 3600) * pi / 180;
        tunnel_slope = 0.01;
    end
end
