classdef BinaryTerrain
% Utility class for writing and reading binary terrain data.

methods (Static)
    function write_from_txt(file_name_no_ext)
        % WRITE_FROM_TXT Write terrain data from the csv file formatted as
        % $station_id, $elevation, $easting, $northing. Units are given in
        % feet, and the order of fields is not changed. Save in single
        % precision if the lidar error is larger than eps.

        fid = fopen([file_name_no_ext '.txt'], 'r');
        data = textscan(fid, '%d %f %f %f', 'HeaderLines', 1, 'Delimiter', ',');
        fclose(fid);

        binary_file = fopen([file_name_no_ext '.bin'], 'w');
        fwrite(binary_file, horzcat(data{2:4}), 'single');
        fclose(binary_file);
    end

    function [X,Y,Elev] = read_file(file_name, total_points, num_points_x, num_points_y)
        % READ_FILE Read in binary data, rearranging fields as necessary
        % and converting to meters. Note that data goes across (left to
        % right) then down, and therefore a transposition is required for
        % reshaping since matlab is column-major.

        assert(num_points_x * num_points_y == total_points);

        fileID = fopen(file_name);
        topo = fread(fileID, [total_points, 3], 'single') * Constants.feet_to_meters;
        fclose(fileID);

        shape = [num_points_x, num_points_y];
        reshaped = map_columns(@(x) reshape(x, shape)', topo, false);
        [Elev, X, Y] = reshaped{:};
    end
end
end

