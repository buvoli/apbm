function [z] = rou(N,rho)
%ROU Returns N rotated roots of unity
% Parameters
%   n   (integer)
%       Number of roots of unity
%   rho (double)
%       Rotation Factor (0<=rho<2*pi/N)
if(nargin==1)
    rho = 0;
end
z = exp(2*pi*1i*(0:N-1)/N + 1i*rho).';
end

