function [amp] = rIMPBDF(z1, z2_vec, options)
%rIMPBDF stability function R(z_1,z_2) = \rho(M(z_1, z_2)) for Polynomial
% IMEX-BDF method.
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   method_handle - (function_handle) handle to desired method
%   options - struct with fields:
%               z - array containing pbm nodes
%               alpha - extrapolation parameter
% RETURNS
%   amp    - amplification factor

% -- set default params ------------------------------------------------------------------------------------------------
if(nargin < 3)
    options = struct();
end

default_value_cell = {
    {'z',     linspace(-1,1,4)}
    {'alpha', 2 / 3}
};
options = setDefaultOptions(options, default_value_cell);

% -- Numerical Parameters ----------------------------------------------------------------------------------------------
z = options.z;
q = length(z);
alpha = options.alpha;

BDFW = transpose(interpW(z, z(end) + alpha, z + alpha));
DDFE = transpose(interpW(z, [], z(end) + alpha));

if( alpha == 2 / (q - 1) && isequal(z, linspace(-1,1,q)) ) % classical IMEX BDF
    BDFW(abs(BDFW) < 1e-14) = 0.0; % remove rounding error from zero entries
end

A = BDFW(1:q, 1:q);
B_im = zeros(q);
B_ex = bsxfun(@times, BDFW(:, q+1), DDFE);
C    = zeros(q);
D_im = [ zeros(q, q-1), BDFW(:, q+1) ];
D_ex = zeros(q);

amp = rIMPBM(z1, z2_vec, alpha, A, B_im, B_ex, C, D_im, D_ex);

end

function w = interpW(x, xp, a)
%QUADW Returns quadrature weights for the nodes x where
%
%   \int^b(i)_a(i) f(x) dx \approx \sum_{j=1}^n w(j,i) f(x(i))
%
% == Parameters ===============================================================
%   x (vector) - nodes
%   a (vector) - left endpoints
%   b (vector) - right endpoints
% =============================================================================

n = length(x);
m = length(xp);
V = zeros(n + m);

for i = 1 : n + m
    for j = 1 : n
        V(i, j) = x(j)^(i-1);
    end
    for j = 1 : m
        V(i, j + n) = (i - 1) * xp(j)^max(0,(i-2));
    end
end

b = (a(:)').^((0:n+m-1)');
w = V \ b;
end