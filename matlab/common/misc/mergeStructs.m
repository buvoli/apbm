function struct_a = mergeStructs(struct_a, struct_b, overwrite)
%SETDEFAULTOPTIONS adds all fields and values of struct_b to struct_a. If overwrite = true, then any duplicated fields 
% in struct_a will be overwritten
% == Parameters ========================================================================================================
%   1. struct_a    (struct) - struct to modify
%   2. struct_b    (struct) - struct that contains additional/updated fields
%   2. overwrite   (bool)   - optional parameter (defaults to true). If overwrite = true, then any duplicated fields 
%                             in struct_a will be overwritten 
% == Returns ===========================================================================================================
%   1. struct_ab   (struct) - struct with updated/empty fields added.
% ======================================================================================================================

if(nargin < 3)
    overwrite = true;
end

fnb = fieldnames(struct_b);
for i = 1 : length(fnb)
    field = fnb{i};
    if(isfield(struct_a, field))
        if(overwrite)
            struct_a.(field) = struct_b.(field);
        end
    else
        struct_a.(field) = struct_b.(field);
    end
end
end