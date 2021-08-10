% ================================================================================
% Produces all stability plots for polynomial BDF5 for the partitioned Dahlquist
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

%% Polynomial BDF Overlay Plots (Colored - Standard Zoom)

Np = 201; % number of points in Re(z) and Im(z)
z1_r    = linspace(0, 6, 101);
z2_real = linspace(-3.5, 1/2, Np);
z2_imag = linspace(-2, 2, Np);

special_contours = @(theta) {
    {6 * exp(1i * theta), {'Color', label_red,    'LineWidth', 2}}
    {3 * exp(1i * theta), {'Color', label_yellow, 'LineWidth', 2}}
    {0 * exp(1i * theta), {'Color',  label_blue,   'LineWidth', 2}}
};
hidden_contours = @(theta) setdiff(exp(1i * theta) * z1_r, [0, 3, 6] * exp(1i * theta));

q0 = 5;
domain = struct( ...
    'theta', [pi, 2 * pi / 3, pi / 2], ...
    'alpha', (2 / (q0 - 1)) * [1, 1/2, 1/4] ...
);

params = struct('q', q0, 'z1_r', z1_r, 'z2_real', z2_real, 'z2_imag', z2_imag, 'symmetric', true, 'hidden_contours', hidden_contours, 'special_contours', special_contours, 'XTicks', [-3 -2 -1 0], 'save_options', save_options);
%multiloopSPMD(@pbdfStabilityPlot, domain, params, 28)  
multiloop(@pbdfStabilityPlot, domain, params)

%% Polynomial BDF Overlay Plots (Colored - Imaginary Axis Zoom)

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

q0 = 5;
domain = struct( ...
    'theta', [2 * pi/3, pi / 2], ...
    'alpha', (2 / (q0 - 1)) * [1, 1/2, 1/4] ...
);

params = struct('q', q0, 'z1_r', z1_r, 'z2_real', z2_real, 'z2_imag', z2_imag, 'symmetric', true, 'hidden_contours', hidden_contours, 'special_contours', special_contours, 'XTicks', [], 'save_options', save_options, 'tag', '-imzoom');
%multiloopSPMD(@pbdfStabilityPlot, domain, params, 28) 
multiloop(@pbdfStabilityPlot, domain, params)

function pbdfStabilityPlot(params)

[t_n, t_d] = rat(params.theta / pi);
[a_n, a_d] = rat(params.alpha);

save_path = fullfile('figures', ['pbdf', num2str(params.q)], ['contours-t-', num2str(t_n), '-', num2str(t_d), '-a-', num2str(a_n), '-', num2str(a_d)]);
if(isfield(params, 'tag'))
    save_path = [save_path, params.tag];
end

amp = @(z1,z2) rIMPBDF(z1, z2, struct('z', linspace(-1,1,params.q), 'alpha', params.alpha));

z1 = exp(1i * params.theta) * params.z1_r;
hidden_contours = params.hidden_contours(params.theta);

fh = OverlayedStabilityPlot(amp, z1, params.z2_real, params.z2_imag, struct('Symmetric', params.symmetric, 'IndexedOverlayEdgeStyle', {params.special_contours(params.theta)}, 'HiddenContours', hidden_contours, 'XTicks', params.XTicks));
exportFigure(fh, params.save_options(save_path))

end