function [results, index_order] = multiloop(f, domain, params)
%SEARCH runs function f for each parameter in the domain.
%   f   (handle)    : f = @(x). Search will locate min or max of f in 
%                     space.
%   domain (struct) : describes search space. of the form
%                     space('var1', linspace(-1,1,10), 'var2', [1,2,5]).
%   params (struct) : OPTIONAL. Can contain additional parameters passed to
%                     f.

% == Set Optional Parameters ============================================ %
if(nargin == 2)
	params = struct();
end
% == Read Space Properties ============================================== %
variables     = fieldnames(domain);
num_variables = length(variables);
dimensions    = zeros(num_variables,1);
for i=1:num_variables
    dimensions(i) = length(domain.(variables{i}));
end

if(nargout >= 1)
    results = cell(prod(dimensions), 1);
end
% == Search Space ======================================================= %
inds   = ones(num_variables, 1);
for i=1:prod(dimensions)
    for j=1:num_variables                                                   % pack variables
        if(iscell(domain.(variables{j})))
        	params.(variables{j}) = domain.(variables{j}){inds(j)};
        else
        	params.(variables{j}) = domain.(variables{j})(inds(j));
        end
    end
    if(nargout >= 1)
        f_eval = f(params);
    else
        f(params);
    end
    if(nargout >= 1)
        results{i} = f_eval;
    end    
    inds = increment(inds, dimensions, num_variables);                      % increment index counter
end
if(nargout >= 1)
    if(length(dimensions)> 1)
        results = reshape(results, dimensions(:)');
    end
    index_order = variables;
end
end

function num = increment(num, bases, num_digits)
%INCREMENT increments "generalized number" num by one.
%   inds        (array) : generalized number. Each digit stored 
%                         incrementally in array.
%   dimensions  (array) : base of each digit position so that
%                         1 <= inds(i) <= bases(i)
%   num_digits  (int)   : integer which is length(num). Added for
%                         performance.
    num(1) = num(1) + 1;
    for i=1:num_digits-1
        if(num(i) > bases(i))
            num(i) = 1;
            num(i+1) = num(i+1) + 1;
        end
    end
end