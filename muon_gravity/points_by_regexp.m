function points = points_by_regexp(table, reg)
    is_included = ~cellfun(@isempty, regexp(table.Properties.RowNames, reg));
    points = table{is_included, Constants.xyz_index}';
end