function visualizeAcousticField(nodes, elements, p, frequencies, conf, mode_idx)
% visualizeAcousticField - Visualize acoustic pressure field
%
% Input:
%   nodes      - Node data structure
%   elements   - Element data structure
%   p          - Pressure field [n_dof × 1] or [n_dof × n_freq]
%   frequencies - Frequency vector [n_freq × 1] (optional)
%   conf       - Configuration structure (optional)
%   mode_idx   - Frequency/mode index to visualize (default: 1)

if nargin < 6
    mode_idx = 1;
end

if nargin < 5
    conf = struct();
    conf.mesher = 'GMSH';
end

if nargin < 4
    frequencies = [];
end

% Determine if p is single field or multiple frequencies
if size(p, 2) > 1
    % Multiple frequencies
    p_plot = p(:, mode_idx);
    if ~isempty(frequencies) && mode_idx <= length(frequencies)
        freq_str = sprintf(' at %.2f Hz', frequencies(mode_idx));
    else
        freq_str = sprintf(' (mode %d)', mode_idx);
    end
else
    % Single field
    p_plot = p;
    if ~isempty(frequencies) && length(frequencies) == 1
        freq_str = sprintf(' at %.2f Hz', frequencies(1));
    else
        freq_str = '';
    end
end

% Get node coordinates
node_coords = nodes.all(:, 1:2);

% Create figure
fig = figure;
set(fig, 'Position', [100, 100, 1000, 800]);

% Extract pressure values at nodes
p_values = zeros(size(node_coords, 1), 1);
for i = 1:size(node_coords, 1)
    dof_idx = nodes.nind(i, 2);
    if dof_idx > 0 && dof_idx <= length(p_plot)
        p_values(i) = real(p_plot(dof_idx));
    end
end

% Plot using contourf
subplot(2, 2, 1);
scatter(node_coords(:, 1), node_coords(:, 2), 50, p_values, 'filled');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Pressure field%s', freq_str));
axis equal;
grid on;

% Plot magnitude
subplot(2, 2, 2);
p_mag = abs(p_plot);
p_mag_values = zeros(size(node_coords, 1), 1);
for i = 1:size(node_coords, 1)
    dof_idx = nodes.nind(i, 2);
    if dof_idx > 0 && dof_idx <= length(p_mag)
        p_mag_values(i) = p_mag(dof_idx);
    end
end
scatter(node_coords(:, 1), node_coords(:, 2), 50, p_mag_values, 'filled');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Pressure magnitude%s', freq_str));
axis equal;
grid on;

% Plot phase
subplot(2, 2, 3);
p_phase = angle(p_plot);
p_phase_values = zeros(size(node_coords, 1), 1);
for i = 1:size(node_coords, 1)
    dof_idx = nodes.nind(i, 2);
    if dof_idx > 0 && dof_idx <= length(p_phase)
        p_phase_values(i) = p_phase(dof_idx);
    end
end
scatter(node_coords(:, 1), node_coords(:, 2), 50, p_phase_values, 'filled');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Pressure phase%s', freq_str));
axis equal;
grid on;

% Plot SPL
subplot(2, 2, 4);
p_ref = 2e-5;
SPL = 20 * log10(abs(p_plot) / p_ref);
SPL_values = zeros(size(node_coords, 1), 1);
for i = 1:size(node_coords, 1)
    dof_idx = nodes.nind(i, 2);
    if dof_idx > 0 && dof_idx <= length(SPL)
        SPL_values(i) = SPL(dof_idx);
    end
end
scatter(node_coords(:, 1), node_coords(:, 2), 50, SPL_values, 'filled');
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Sound Pressure Level (SPL)%s', freq_str));
axis equal;
grid on;

sgtitle(sprintf('Acoustic Field Visualization%s', freq_str), 'FontSize', 14);

end

