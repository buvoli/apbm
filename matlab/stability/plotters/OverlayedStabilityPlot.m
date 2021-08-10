function [fh, data_raw] = OverlayedStabilityPlot(amp, z1, z2_real, z2_imag, options)
%OVERLAYEDSTABILITYDATA produces the raw data for an overlayed stability plot
%   amp     (handle) - function of two arguments @(z1 - scalar, z2 - vector) producing amp factors
%   z1      (scalar) - z1 = lambda_1 * h
%   z2_real (vector) - real part of z2 = lambda_2 * h. 
%   z2_imag (vector) - imaginary part of z2 = z2 = lambda_2 * h

if(nargin == 4)
    options = struct();
end

data_raw = OverlayedStabilityData(amp, z1, z2_real, z2_imag, options);
fh       = OverlayedStabilityPlotter(data_raw, z1, z2_real, z2_imag, options);
end