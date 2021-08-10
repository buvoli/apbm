% ================================================================================
% Produces all stability plots for polynomial Radau4 for the partitioned Dahlquist
% ================================================================================

% -- include library files -------------------------------------------------------
addpath(genpath('../../stability'));
addpath(genpath('../../common'));

%% === Set Parameters ===
clear; clc; close all;

% Global Plot Parameters
label_red    = [214 7 37] / 255;
label_yellow = [253, 209, 7] / 255;
label_blue   = [53 146 177] / 255;
label_black  = [.1 .1 .1];

save_options = @(save_path) struct('SavePath', save_path, 'PaperPosition', [0 0 8 7], 'Format', 'pdf');

%% Polynomial Fully Implicit Radau (Regular Zoom)

Np = 201; Lz = 2.2; % number of points in Re(z) and Im(z)
z1_r    = linspace(0, 6, 101);
z2_real = linspace(-2.2, 2, Np);
z2_imag = linspace(-2, 2, Np);

special_contours = @(theta) {
    {6 * exp(1i * theta), {'Color', label_red,    'LineWidth', 2}}
    {3 * exp(1i * theta), {'Color', label_yellow, 'LineWidth', 2}}
    {0 * exp(1i * theta), {'Color',  label_blue,   'LineWidth', 2}}
};
hidden_contours = @(theta) setdiff(exp(1i * theta) * z1_r, [0, 3, 6] * exp(1i * theta));

q0 = 4;
domain = struct( ...
    'theta', [pi, 2 * pi / 3, pi / 2], ...
    'kappa', [0 1 2 3 4] ...
);

params = struct('q', q0, 'alpha', 2, 'z1_r', z1_r, 'z2_real', z2_real, 'z2_imag', z2_imag, 'symmetric', true, 'special_contours', special_contours, 'hidden_contours', hidden_contours, 'save_options', save_options, 'XTicks', [-2 -1 0 1 2]);
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', false)))
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', true)))

%% Polynomial Fully Implicit Radau (Magnified)

Np = 201; Lz = 0.5; % number of points in Re(z) and Im(z)
z1_r    = linspace(0, 6, 101);
z2_real = linspace(-Lz, Lz, Np);
z2_imag = linspace(-Lz, Lz, Np);

special_contours = @(theta) {
    {6 * exp(1i * theta), {'Color', label_red,    'LineWidth', 2}}
    {3 * exp(1i * theta), {'Color', label_yellow, 'LineWidth', 2}}
    {0 * exp(1i * theta), {'Color',  label_blue,   'LineWidth', 2}}
};
hidden_contours = @(theta) setdiff(exp(1i * theta) * z1_r, [0, 3, 6] * exp(1i * theta));

q0 = 4;
domain = struct( ...
    'theta', [pi, 2 * pi / 3, pi / 2], ...
    'kappa', [0 1 2 3 4] ...
);

params = struct('q', q0, 'alpha', 2, 'z1_r', z1_r, 'z2_real', z2_real, 'z2_imag', z2_imag, 'symmetric', true, 'special_contours', special_contours, 'hidden_contours', hidden_contours, 'save_options', save_options, 'XTicks', [], 'tag', '-zoom');
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', false)))
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', true)))

%% Polynomial Fully Implicit Radau (Imaginary Axis Zoom)

Np = 201; Lz = .1; % number of points in Re(z) and Im(z)
z1_r    = linspace(0, 6, 101);
z2_real = linspace(-2*Lz + Lz/4, Lz/4, Np);
z2_imag = linspace(-Lz, Lz, Np);

special_contours = @(theta) {
    {6 * exp(1i * theta), {'Color', label_red,    'LineWidth', 2}}
    {3 * exp(1i * theta), {'Color', label_yellow, 'LineWidth', 2}}
    {0 * exp(1i * theta), {'Color',  label_blue,   'LineWidth', 2}}
};
hidden_contours = @(theta) setdiff(exp(1i * theta) * z1_r, [0, 3, 6] * exp(1i * theta));

q0 = 4;
domain = struct( ...
    'theta', [2 * pi / 3, pi / 2], ...
    'kappa', [0 1 2 3 4] ...
);

params = struct('q', q0, 'alpha', 2, 'z1_r', z1_r, 'z2_real', z2_real, 'z2_imag', z2_imag, 'symmetric', true, 'special_contours', special_contours, 'hidden_contours', hidden_contours, 'save_options', save_options, 'XTicks', [], 'tag', '-imzoom');
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', false)))
multiloop(@pRadauStabilityPlot, domain, mergeStructs(params, struct('star', true)))

function pRadauStabilityPlot(params)

[t_n, t_d] = rat(params.theta / pi);
[a_n, a_d] = rat(params.alpha);

if(params.star)
    amp = radauFIS(params.q, params.alpha, params.kappa);
    method_name = 'pRadauS';
else
    amp = radauFI(params.q, params.alpha, params.kappa);
     method_name = 'pRadau';
end

save_path = fullfile('figures', [method_name, num2str(params.q)], ['contours-t-', num2str(t_n), '-', num2str(t_d), '-a-', num2str(a_n), '-', num2str(a_d), '-k-', num2str(params.kappa)]);
if(isfield(params, 'tag'))
    save_path = [save_path, params.tag];
end

z1 = exp(1i * params.theta) * params.z1_r;
hidden_contours = params.hidden_contours(params.theta);

fh = OverlayedStabilityPlot(amp, z1, params.z2_real, params.z2_imag, struct('Symmetric', params.symmetric, 'IndexedOverlayEdgeStyle', {params.special_contours(params.theta)}, 'HiddenContours', hidden_contours, 'XTicks', params.XTicks));
exportFigure(fh, params.save_options(save_path))

end

function amp = radauFIS(q, alpha, kappa)
    amp = @(z1,z2) rFIP(z1, z2, struct('z', [-1; -flip(radaupts(q-1))], 'alpha', alpha, 'b_ind', q * ones(q, 1), 'I_AODS', 2 : q, 'E_AIDS', 1 : q, 'kappa', kappa));
end

function amp = radauFI(q, alpha, kappa)
    amp = @(z1,z2) rFIP(z1, z2, struct('z', [-1; -flip(radaupts(q-1))], 'alpha', alpha, 'b_ind', q * ones(q, 1), 'I_AODS', 2 : q, 'E_AIDS', 2 : q, 'kappa', kappa));
end