function [x, y] = approximatePositionVector(xi, eta, nodeCoords, elementType)
% approximatePositionVector - Approximate continuous position vector using shape functions
%
% Input:
%   xi          - Natural coordinate xi (scalar or vector)
%   eta         - Natural coordinate eta (scalar or vector, for 2D elements only)
%   nodeCoords  - Node coordinates matrix [n_nodes × 2] or [n_nodes × 3]
%                 Each row: physical coordinates (x, y) or (x, y, z)
%                 Node order must follow GMSH convention
%   elementType - Element type: 1 (Line2), 8 (Line3), 3 (Quadrangle4), 10 (Quadrangle9)
%
% Output:
%   x  - Physical coordinate x [n_points × 1]
%   y  - Physical coordinate y [n_points × 1]
%
% Isoparametric mapping:
%   x(xi, eta) = sum of N_i(xi, eta) · x_i
%   y(xi, eta) = sum of N_i(xi, eta) · y_i
%   where N_i are shape functions, x_i and y_i are physical coordinates of node i

% Check inputs
if nargin < 4
    error('Need 4 input arguments: xi, eta, nodeCoords, elementType');
end

% Convert to column vectors
xi = xi(:);
if nargin >= 2 && ~isempty(eta)
    eta = eta(:);
end

% Check node coordinates dimension
if size(nodeCoords, 2) < 2
    error('nodeCoords must contain at least x and y coordinates');
end

% Compute shape functions
if elementType == 1 || elementType == 8
    % 1D element
    [N, ~] = shapeFunctions(xi, [], elementType);
else
    % 2D element
    if nargin < 2 || isempty(eta)
        error('2D elements require eta parameter');
    end
    [N, ~] = shapeFunctions(xi, eta, elementType);
end

% Check node count match
n_nodes_expected = size(N, 2);
n_nodes_provided = size(nodeCoords, 1);

if n_nodes_expected ~= n_nodes_provided
    error('Node count mismatch: expected %d nodes, provided %d nodes', ...
        n_nodes_expected, n_nodes_provided);
end

% Compute continuous position vector
% x(xi, eta) = sum of N_i(xi, eta) · x_i
% y(xi, eta) = sum of N_i(xi, eta) · y_i
x = N * nodeCoords(:, 1);
y = N * nodeCoords(:, 2);

end
