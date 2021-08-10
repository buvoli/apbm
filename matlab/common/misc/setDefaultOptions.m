function options = setDefaultOptions(options, default_field_value_pairs)
%SETDEFAULTOPTIONS fills in empty fields in struct with defailt values
% == Parameters ========================================================================================================
%   1. options                      (struct) - struct to modify
%   2. default_field_value_pairs    (cell)   - 2xn cell where default_field_value_pairs{i} = {key, default_value}
% == Returns ===========================================================================================================
%   1. options                      (struct) - struct with empty fields added.
% ======================================================================================================================

for i = 1 : length(default_field_value_pairs)
    field = default_field_value_pairs{i}{1};
    if(~isfield(options, field))
        value = default_field_value_pairs{i}{2};
        options.(field) = value;
    end
end

end