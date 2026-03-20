function [p, frequencies] = computeFrequencyResponse(K, M, f, freq_range, air, varargin)
% computeFrequencyResponse - Compute frequency response
%
% Input:
%   K          - Stiffness matrix [n_dof × n_dof] (with BC applied)
%   M          - Mass matrix [n_dof × n_dof] (with BC applied)
%   f          - Load vector [n_dof × 1] OR source info structure
%   freq_range - Frequency range [f_min, f_max, n_points] or frequency vector
%   air        - Air parameters structure (contains speed of sound c)
%   varargin   - Optional: nodes structure, boundary_nodes, boundary_type
%                If provided, boundary conditions will be applied to load vector
%
%   If f is a structure with fields 'source_position', 'source_strength',
%   the load vector will be recomputed for each frequency.
%
% Output:
%   p          - Pressure solution [n_dof × n_freq]
%   frequencies - Frequency vector [n_freq × 1] (Hz)

if nargin < 5 || ~isfield(air, 'c')
    error('Need air structure with speed of sound c');
end

if nargin < 4
    freq_range = [10, 200, 100];
end

% Parse frequency range
if length(freq_range) == 3
    f_min = freq_range(1);
    f_max = freq_range(2);
    n_points = freq_range(3);
    frequencies = linspace(f_min, f_max, n_points);
else
    frequencies = freq_range(:);
    n_points = length(frequencies);
end

n_dof = size(K, 1);
p = zeros(n_dof, n_points);

% Check if f is a structure (source info) or a vector
use_frequency_dependent_source = isstruct(f) && isfield(f, 'source_position');

% Check for boundary condition info
apply_bc_to_load = false;
if nargin >= 6
    nodes = varargin{1};
    if nargin >= 7 && ~isempty(varargin{2})
        boundary_nodes = varargin{2};
        if nargin >= 8 && ~isempty(varargin{3})
            boundary_type = varargin{3};
        else
            boundary_type = 'hard';
        end
        apply_bc_to_load = true;
    end
end

fprintf('Computing frequency response for %d frequencies...\n', n_points);

for i = 1:n_points
    freq = frequencies(i);
    omega = 2 * pi * freq;
    k = omega / air.c;
    
    % Recompute load vector if source info is provided
    if use_frequency_dependent_source
        if ~exist('nodes', 'var')
            error('Need nodes structure when using source info');
        end
        f_current = addAcousticSource(nodes, f.source_position, ...
            f.source_strength, omega, air);
    else
        f_current = f;
    end
    
    % Apply boundary conditions to load vector if needed
    if apply_bc_to_load
        [~, ~, f_current] = applyBoundaryConditions(K, M, f_current, ...
            boundary_nodes, boundary_type, nodes);
    end
    
    % Build system matrix: A = K - k²M
    A = K - (k^2) * M;
    
    % Solve linear system: A·p = f
    try
        p(:, i) = A \ f_current;
    catch ME
        warning('Failed to solve at frequency %.2f Hz: %s', freq, ME.message);
        p(:, i) = NaN(n_dof, 1);
    end
    
    if mod(i, max(1, floor(n_points/10))) == 0
        fprintf('  Progress: %d/%d frequencies\n', i, n_points);
    end
end

fprintf('Done\n');

end

