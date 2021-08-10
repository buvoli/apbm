function amp = rIMCPBM(z1, z2_vec, methods)
%RIMPBM Stability function for composite implicit-explicit polynomial PBM.
%   
%   y^[n+1] = M_n \circ M_{n-1} ... M_1 y^[n]
%
%   where each M_j is a PBM of the form:
%
%   y^[n+1] = A * y^[n] + B_im * f_im(y^[n]) + B_ex * f_ex(y^[n]) + C *
%             + C * y^[n+1] + D_im * f_im(y^[n+1]) + D_ex * f_ex(y^[n+1])
%
% PARAMETERS
%   z1     (scalar) - implicit term f_im(y) = z_1 * y where z_1 = h * \lambda_1.
%   z2     (vector) - explicit term f_ex(y) = z_2 * y where z_2 = h * \lambda_2. Multiple z_2 can be passed in as a vector
%   methods (struct) - struct array containing fields: A, B_im, B_ex, C,
%   D_im, D_ex, kappa

num_methods = length(methods);
q = size(methods(1).A,1);

num_z2 = length(z2_vec);
amp    = zeros(num_z2, 1);
for i = 1 : num_z2    
    z2 = z2_vec(i);
    M  = eye(q);
    for j = 1 : num_methods
        M = ((eye(q) - methods(j).C - z1 * methods(j).D_im - z2 * methods(j).D_ex) \ (methods(j).A + z2 * methods(j).B_ex + z1 * methods(j).B_im))^methods(j).kappa * M;
    end
    amp(i) = max(abs(eig(M)));
end

end