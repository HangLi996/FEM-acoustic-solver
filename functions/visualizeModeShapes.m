function visualizeModeShapes(nodes, elements, eigenvectors, mode_idx, eigenfrequencies, conf)
% visualizeModeShapes - Visualize pressure mode shapes
%
% Input:
%   nodes            - Node data structure
%   elements         - Element data structure
%   eigenvectors     - Eigenvector matrix [n_dof × n_modes]
%   mode_idx         - Mode indices to visualize (can be vector for multiple modes)
%   eigenfrequencies - Eigenfrequencies [n_modes × 1] (Hz)
%   conf             - Configuration structure (optional, for plotMesh)
%
% Output:
%   None, displays figure window

if nargin < 6
    conf = struct();
    conf.mesher = 'GMSH';
end

if nargin < 5 || isempty(eigenfrequencies)
    eigenfrequencies = NaN(size(eigenvectors, 2), 1);
end

if isscalar(mode_idx)
    mode_idx = mode_idx;
end

n_modes_available = size(eigenvectors, 2);
mode_idx = mode_idx(mode_idx <= n_modes_available);
mode_idx = mode_idx(mode_idx > 0);

if isempty(mode_idx)
    error('Invalid mode indices');
end

fig = figure;
set(fig, 'Position', [100, 100, 1200, 800]);

n_plots = length(mode_idx);
n_cols = min(3, n_plots);
n_rows = ceil(n_plots / n_cols);

for plot_idx = 1:n_plots
    mode = mode_idx(plot_idx);
    
    mode_shape = eigenvectors(:, mode);
    
    if ~isnan(eigenfrequencies(mode))
        freq = eigenfrequencies(mode);
        title_str = sprintf('Mode %d: f = %.2f Hz', mode, freq);
    else
        title_str = sprintf('Mode %d', mode);
    end
    
    subplot(n_rows, n_cols, plot_idx);
    hold on;
    
    try
        plotMesh(gca, nodes, elements, conf, [], [1 1 1], 'k', [0 90], 'northwest');
    catch
    end
    
    node_coords = nodes.all(:, 1:2);
    
    scatter(node_coords(:, 1), node_coords(:, 2), 50, mode_shape, 'filled');
    
    colorbar;
    colormap('jet');
    
    title(title_str, 'FontSize', 12);
    xlabel('x (m)', 'FontSize', 10);
    ylabel('y (m)', 'FontSize', 10);
    axis equal;
    grid on;
end

sgtitle('Pressure mode shapes', 'FontSize', 14, 'FontWeight', 'bold');

end

