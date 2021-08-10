function [results, index_order] = multiloopSPMD(f, domain, params, num_procs, partitioning)
%SEARCH runs function f for each parameter in the domain.
%   f   (handle)    : f = @(x). Search will locate min or max of f in 
%                     space.
%   domain (struct) : describes search space. of the form
%                     space('var1', linspace(-1,1,10), 'var2', [1,2,5]).
%   params (struct) : OPTIONAL. Can contain additional parameters passed to
%                     f.

% == Set Optional Parameters ============================================ %
if(nargin <= 2)
	params = struct();
end
if(nargin <= 3)
    num_procs = 1;
end
if(nargin <= 4)
    partitioning = 'interleaved';
end
% == Read Space Properties ============================================== %
variables     = fieldnames(domain);
num_variables = length(variables);
dimensions    = zeros(num_variables,1);
for i=1:num_variables
    dimensions(i) = length(domain.(variables{i}));
end

% == Search Space ======================================================= %
num_tasks = prod(dimensions);
store_output = nargout >= 1;
num_procs = min(num_tasks, num_procs);
spmd(num_procs)
    
    indices = parTaskSplit(num_tasks, num_procs, labindex, partitioning);
    results = cell(1, length(indices));
    count = 1;
    
    for i = indices
        gnum = ind2Gnum(i, dimensions, num_variables, num_tasks);
        params.labindex = labindex;
        for j=1:num_variables                                                 % pack variables
            if(iscell(domain.(variables{j})))
                params.(variables{j}) = domain.(variables{j}){gnum(j)};
            else
                params.(variables{j}) = domain.(variables{j})(gnum(j));
            end
        end
        if(store_output)
            results{count} = f(params);
            count = count + 1;
        else
            f(params);
        end
    end  
end
results = [ results{:} ]; % merge all results into one cell array
if(nargout >= 1)
    results = results(reshuffleInds(num_tasks, num_procs, partitioning));    
    if(length(dimensions)> 1)
        results = reshape(results, dimensions(:)');
    end
    index_order = variables;
end
end

function gnum = ind2Gnum(index, bases, num_digits, b_prod)
%indToGnum converts and index to a generalized number
% index         (integer)   : index of number
% bases         (vector)    : bases of each generalized digit
% num_digits    (integer)   : should be equal to length(bases)
% b_prod        (integer)   : should be equal to prod(bases)

index = index - 1;  % shift by one since 0 is really the first index
gnum = zeros(1, num_digits);
for i = 1 : num_digits
    b_prod = b_prod / bases(i); 
    gnum(i) = floor( index / b_prod );
    index = index - gnum(i) * b_prod;
end
gnum = gnum + 1; % shift by one since Matlab indices start at one

end

function indices = parTaskSplit(num_tasks, num_proc, proc_index, partitioning)
%PARTASKSPLIT Summary of this function goes here
%   num_proc (index 1 .... to numproc)
%   num_tasks (indexed 1 ... to num_tasks)
%   patitioning (str) 'contiguous' or 'block'
    
    if(nargin < 4)
        partitioning = 'contiguous';
    end

    if((nargin == 2) || isempty(proc_index))
        indices = arrayfun(@(pi) parTaskSplit(num_tasks, num_proc, pi, partitioning), 1:num_proc, 'UniformOutput', false);
    else
        if(strcmp(partitioning, 'contiguous'))       
            rem = mod(num_tasks, num_proc);
            dlt = floor(num_tasks/num_proc);
            if((proc_index - 1) < rem)
                start_correction = proc_index;
                stop_correction  = 1;
            else
                start_correction = rem + 1;
                stop_correction  = 0;
            end
            start = dlt * (proc_index - 1) + start_correction;
            stop  = start + (dlt - 1) + stop_correction;
            indices = start:stop;
        elseif(strcmp(partitioning, 'interleaved'))
            indices = proc_index:num_proc:num_tasks;
        end
    end
end

function inds = reshuffleInds(num_tasks, num_procs, partitioning) 
    indices = [];
    for labindex = 1 : num_procs
        indices = [indices, parTaskSplit(num_tasks, num_procs, labindex, partitioning)];
    end
    [~, inds] = sort(indices);
end