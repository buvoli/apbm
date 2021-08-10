function [data_raw] = OverlayedStabilityData(amp, z1, z2_real, z2_imag, options)
%OVERLAYEDSTABILITYDATA produces the raw data for an overlayed stability plot
%   amp     (handle) - function of two arguments @(z1 - scalar, z2 - vector) producing amp factors
%   z1      (scalar) - z1 = lambda_1 * h
%   z2_real (vector) - real part of z2 = lambda_2 * h. 
%   z2_imag (vector) - imaginary part of z2 = z2 = lambda_2 * h
%   symmetric (bool) - if set to true, stability function is max(R(z1,z2), R(z1^*,z2))

default_options = {
    {'Symmetric', false}
};
options = setDefaultOptions(options, default_options);

num_z1   = min(length(z1));
data_raw = cell(length(num_z1));
[X, Y]   = meshgrid(z2_real, z2_imag);
sG       = size(X);
z2       = X(:) + 1i * Y(:);

for i = 1 : num_z1
    if(options.Symmetric && imag(z1(i)) ~= 0)
       a = max(amp(z1(i), z2), amp(conj(z1(i)), z2)); 
    else
        a = amp(z1(i), z2);
    end
    data_raw{i} = reshape(a, sG);
end
end