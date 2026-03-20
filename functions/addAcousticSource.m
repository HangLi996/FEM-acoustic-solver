function f = addAcousticSource(nodes, source_position, source_strength, omega, air)
% addAcousticSource - Add point acoustic source to load vector
%
% Input:
%   nodes          - Node data structure
%   source_position - Source position [x, y] or [x, y, z]
%   source_strength - Source strength Q (volume velocity)
%   omega          - Angular frequency [rad/s]
%   air            - Air parameters structure (contains rho)
%
% Output:
%   f - Load vector [n_dof × 1]

if nargin < 5 || ~isfield(air, 'rho')
    error('Need air structure with density rho');
end

if nargin < 4
    error('Need angular frequency omega');
end

if nargin < 3
    source_strength = 1.0;
end

% Find closest node to source position
node_coords = nodes.all(:, 1:2);
n_nodes = size(node_coords, 1);
distances = sqrt(sum((node_coords - repmat(source_position(1:2), n_nodes, 1)).^2, 2));
[~, closest_node_idx] = min(distances);

% Get DOF index
source_dof = nodes.nind(closest_node_idx, 2);

% Get actual DOF count from nodes structure
% Find max DOF index to determine size
max_dof = max(nodes.nind(:, 2));
n_dof = max_dof;

% Initialize load vector with correct size
f = sparse(n_dof, 1);

if source_dof > 0 && source_dof <= n_dof
    % Add point source: f = -j*rho*omega*Q
    f(source_dof) = -1j * air.rho * omega * source_strength;
else
    warning('Invalid DOF index %d for source node %d (max DOF: %d)', ...
        source_dof, closest_node_idx, n_dof);
end

end

