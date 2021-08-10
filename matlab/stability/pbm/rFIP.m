function [amp] = rFIP(z1, z2_vec, options)
%rFIP stability function R(z_1,z_2) = \rho(M(z_1, z_2)) for Polynomial IMEX fully implicit method.
% PARAMETERS
%   z1     (scalar) - exponential term: z_1 = h * \lambda_1
%   z2     (vector) - exponential term: z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   method_handle - (function_handle) handle to desired method
%   options - struct with fields:
%               z      -  array containing pbm nodes
%               b_ind  - indices for left endpoints
%               E_AIDS - active input derivative set for explicit nonlinear term
%               I_AODS - active input derivative set for linear term
%               alpha  - extrapolation parameter
%               kappa  - number of correction interations
% RETURNS
%   amp    - amplification factor

% -- set default params ---------------------------------------------------
if(nargin < 3)
    options = struct();
end

default_value_cell = {
    {'z',      linspace(-1,1,3)}
    {'b_ind',  ones(1,3)}
    {'E_AIDS', 1:3}
    {'I_AODS', 1:3}
    {'alpha',  2}
    {'kappa',  10}
};
options = setDefaultOptions(options, default_value_cell);

% -- Numerical Parameters -------------------------------------------------
z = options.z;
q = length(z);
alpha = options.alpha;
b = z(options.b_ind);   % left endpoints
kappa = options.kappa;  % iterations

% -- Implicit Weights -----------------------------------------------------
I_AODS = options.I_AODS; % active output derivative indices for linear term

% ----> propagator
AP = zeros(q);
AP(1:q, options.b_ind) = 1;
BIP = zeros(q);
BIP(:, I_AODS) = transpose(quadW(z(I_AODS), b - alpha, z));

% ----> interator (always use left endpoint)
AI = zeros(q);
AI(1:q, 1) = 1;
BII = zeros(q);
BII(:, I_AODS) = transpose(quadW(z(I_AODS), z(1), z));

% -- Explicit weights -----------------------------------------------------
E_AIDS = options.E_AIDS; % active input derivative indices for linear term

% ----> propagator
BEP = zeros(q);
BEP(:, E_AIDS) = transpose(quadW(z(E_AIDS), b, z + alpha));
% ----> interator (reuse implicit iterator)
BEI = BII;

ZM = zeros(q);
methods = struct('A', {AP, AI}, 'B_im', {ZM, ZM}, 'B_ex', {BEP, BEI}, 'C', {ZM, ZM}, 'D_im', {BIP, BII}, 'D_ex', {ZM,ZM}, 'kappa', {1, kappa});
amp = rIMCPBM(z1, z2_vec, methods);

end

function w = quadW(x, a, b)
%QUADW Returns quadrature weights for the nodes x where
%
%   \int^b(i)_a(i) f(x) dx \approx \sum_{j=1}^n w(j,i) f(x(i))
%
% == Parameters ===============================================================
%   x (vector) - nodes
%   a (vector) - left endpoints
%   b (vector) - right endpoints
% =============================================================================

V = transpose(fliplr(vander(x)));
p = (1:length(x))';
b = ((b(:)').^p - (a(:)').^p ) ./ p;
w = V \ b;
end