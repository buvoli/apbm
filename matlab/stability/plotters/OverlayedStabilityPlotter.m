function [figure_handle] = OverlayedStabilityPlotter(data_raw, z1, z2_real, z2_imag, options)
%OVERLAYEDSTABILITYDATA produces the raw data for an overlayed stability plot
%   data_raw (struct) - data calculated using OverlayedStabilityData
%   z1       (scalar) - z1 = lambda_1 * h
%   z2_real  (vector) - real part of z2 = lambda_2 * h. 
%   z2_imag  (vector) - imaginary part of z2 = z2 = lambda_2 * h

if(nargin == 4)
    options = struct();
end
default_options = {
    {'FigureIndex',             1}
    {'ClearFigure',             true}
    {'DrawAxis',                true}
    {'Overlay',                 'all'}
    {'OverlayEdgeColor',        .3 * [1 1 1]}
    {'OverlayLineWidth',        1}
    {'Union',                   'all'}
    {'UnionEdgeColor',          .3 * [1 1 1]}
    {'UnionBackgroundColor',    .8 * [1 1 1]} 
    {'UnionLineWidth',          1}
    {'IndexedOverlayEdgeStyle', {}}  % nx2 cell array A{j,:} = {z1_special, line_style}. contour with z1=z1_special will be created using arguments line_style{:}
    {'HiddenContours',          []}  % contours for any z1 included in this array will not be shown (but will be used to compute Union
    {'ContourLevels'            [0 1]}
    {'FontSize',                12}
    {'FontName',                'Minion Pro'}
    {'AxisSquare',              true}
    {'NumTicks',                5}
    {'XTicks',                  []} % can be used to ovverride NumTicks
    {'YTicks',                  []} % can be used to ovverride NumTicks
};
options = setDefaultOptions(options, default_options);


% -- Set Figure Index --------------------------------------------------------------------------------------------------
if(isempty(options.FigureIndex))
    figure_handle = figure();
else
    figure_handle = figure(options.FigureIndex);
    if(options.ClearFigure)
        clf;
    end    
end


% -- create overlay plot -----------------------------------------------------------------------------------------------
if(strcmp(options.Union, 'all'))
    union_inds = 1 : length(data_raw);
else
    union_inds = options.Overlay;
end
num_union = length(union_inds);

if(num_union > 0)
    union = data_raw{1};
    for i = 2 : num_union
        union = max(union, data_raw{i});
    end
    if(min(union(:)) < max(options.ContourLevels)) % only plot if at least one value is below threshold
        contourf(z2_real, z2_imag, union, options.ContourLevels, 'color', options.UnionEdgeColor, 'linewidth', options.UnionLineWidth);
        hold on;
        if(~isempty(options.UnionBackgroundColor))
            colormap([options.UnionBackgroundColor; 1 1 1]);
        end
    else
        warning('there were not stable points in union set');
    end
end

% -- draw axis ---------------------------------------------------------------------------------------------------------
if(options.DrawAxis)
    plot([min(z2_real) max(z2_real)], [0 0], 'k'); hold on; 
    plot([0 0], [min(z2_imag) max(z2_imag)], 'k');
end

% -- create overlay plot -----------------------------------------------------------------------------------------------
num_contours = length(data_raw);
if(isempty(options.IndexedOverlayEdgeStyle))
    special_inds = [];
else
    special_inds = cellfun(@(c) c{1}, options.IndexedOverlayEdgeStyle); % special inds that need to be colored
end
hidden_inds = options.HiddenContours;

[~, inds] = setdiff(z1, hidden_inds); % only show non-hidden inds
inds = transpose(inds(:)); % turn into row vector

for i = inds    
    [special_color, ind] = ismember(z1(i), special_inds); % -- see if linestyle is special, otherwise use default
    if(special_color)
        line_style = options.IndexedOverlayEdgeStyle{ind}{2};
    else
        line_style = {'color', options.OverlayEdgeColor, 'linewidth', options.OverlayLineWidth};
    end
    contour(z2_real, z2_imag, data_raw{i}, options.ContourLevels, line_style{:}); hold on;
end

% -- set Font and axis ----
set(gca, 'FontSize', options.FontSize, 'FontName', options.FontName);

if(options.AxisSquare)
    axis square;
end

if(isempty(options.XTicks))
    xtick_vec = linspace(min(z2_real), max(z2_real), options.NumTicks);
    xticks(xtick_vec);
else
    xticks(options.XTicks);
end

if(isempty(options.YTicks))
    ytick_vec = linspace(min(z2_imag), max(z2_imag), options.NumTicks);
    yticks(ytick_vec);
else
    yticks(options.YTicks);
end

hold off;
end