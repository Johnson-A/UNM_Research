function [table, measured_points] = build_table()
    table = readtable('point_database.csv', 'ReadRowNames', true);
    table = sortrows(table, 'RowNames');
    table.Properties.VariableUnits = {'', '', '', ''};
    
    % TODO: Comment NAD83 State plane coordinates
    
    stations = textscan(fopen('gz_measurements.csv', 'r'), '%s %f64 %f64', 'Delimiter',',');
    
    all_names = stations{1};
    all_errors = stations{2};
    all_measurements = stations{3};
    measured_points = unique(all_names);
    
    measurement_table = cell2table(cell(size(table, 1), 2), ...
        'VariableNames', {'Measurements', 'Errors'}, ....
        'RowNames', table.Properties.RowNames);
    
    for name = measured_points',
        is_a_measurement = strcmp(all_names, name);
        assert(any(strcmp(table.Properties.RowNames, name)), ['Station ' name ' is not in the database of points'])
        % TODO: Is this the write way to assign (indexing gives 1x1 cell)
        measurement_table{name, 1:2} = {all_measurements(is_a_measurement), all_errors(is_a_measurement)};
    end
    
    table = [table, measurement_table];
end