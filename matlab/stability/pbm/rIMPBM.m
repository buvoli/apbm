function amp = rIMPBM(z1, z2_vec, alpha, A, B_im, B_ex, C, D_im, D_ex)
%rIMPBM returns stability function magnitude for implict-explicit
%polynomial PBM
%   
%   y^[n+1] = A * y^[n] + B_im * f_im(y^[n]) + B_ex * f_ex(y^[n]) + C *
%             + C * y^[n+1] + D_im * f_im(y^[n+1]) + D_ex * f_ex(y^[n+1])
%
% PARAMETERS
%   z1     (scalar) - implicit term f_im(y) = z_1 * y where z_1 = h * \lambda_1.
%   z2     (vector) - explicit term f_ex(y) = z_2 * y where z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector

q = size(A,1);

num_z2 = length(z2_vec);
amp    = zeros(num_z2, 1);
z1     = z1 / alpha; % z1 -> r * \lambda_1
for i = 1 : num_z2    
    z2 = z2_vec(i) / alpha;
    M = (eye(q) - C - z1 * D_im - z2 * D_ex) \ (A + z2 * B_ex + z1 * B_im);
    amp(i) = max(abs(eig(M)));
end

end