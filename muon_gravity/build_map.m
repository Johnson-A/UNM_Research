function [measurement_map, uniq_str_ids, locations] = build_map()
    stations = textscan(fopen('gz_measurements.csv','r'), '%s %f64 %f64', 'Delimiter',',');
    
    all_names = stations{1};
    all_errors = stations{2};
    all_values = stations{3};
    
    uniq_str_ids = sort(unique(all_names));
    measurement_map = containers.Map();

    num_pts = numel(uniq_str_ids);
    locations = zeros(3, num_pts);
    
    for pt = 1:num_pts,
        name = cell2mat(uniq_str_ids(pt));
        is_a_measurement = strcmp(all_names, name);
        measurement_map(name) = [all_errors(is_a_measurement), all_values(is_a_measurement)];
        
        in_pt_database = strcmp(Constants.pt_names, name);
        assert(any(in_pt_database), ['Error: The point ' name ' is not in the database']);
        assert(sum(in_pt_database) == 1, ['Error: There are multiple points named ' name ' in the database']);
        locations(:, pt) = Constants.all_pts(:, strcmp(Constants.pt_names, name));
    end
end