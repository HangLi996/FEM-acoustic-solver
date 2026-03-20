function [J, detJ] = computeJacobian(xi, eta, nodeCoords, elementType)
% computeJacobian - Compute Jacobian matrix and determinant at natural coordinates
%
% Input:
%   xi, eta      - Natural coordinates (scalar or vector, typically Gauss points)
%   nodeCoords  - Element node coordinates [n_nodes × 2]
%   elementType - Element type: 1 (Line2), 8 (Line3), 3 (Quadrangle4), 10 (Quadrangle9)
%
% Output:
%   J    - Jacobian matrix [2 × 2] for 2D, [1 × 1] for 1D
%   detJ - Jacobian determinant [scalar]

if nargin < 4
    error('Need 4 input arguments: xi, eta, nodeCoords, elementType');
end

xi = xi(:);
if nargin >= 2 && ~isempty(eta)
    eta = eta(:);
end

is1D = (elementType == 1 || elementType == 8);
is2D = (elementType == 3 || elementType == 10);

if ~is1D && ~is2D
    error('Unsupported element type: %d', elementType);
end

if is2D && (nargin < 2 || isempty(eta))
    error('2D elements require eta parameter');
end

if is1D
    [N, dN] = shapeFunctions(xi, [], elementType);
else
    [N, dN] = shapeFunctions(xi, eta, elementType);
end

n_nodes_expected = size(N, 2);
n_nodes_provided = size(nodeCoords, 1);

if n_nodes_expected ~= n_nodes_provided
    error('Node count mismatch: expected %d, provided %d', ...
        n_nodes_expected, n_nodes_provided);
end

n_points = length(xi);
if is1D
    J = zeros(1, 1, n_points);
    detJ = zeros(n_points, 1);
else
    J = zeros(2, 2, n_points);
    detJ = zeros(n_points, 1);
end

for p = 1:n_points
    if is1D
        dN_dxi = squeeze(dN(p, :, 1));
        if size(dN_dxi, 1) > size(dN_dxi, 2)
            dN_dxi = dN_dxi';
        end
        
        J(1, 1, p) = dN_dxi * nodeCoords(:, 1);
        detJ(p) = abs(J(1, 1, p));
    else
        dN_dxi = squeeze(dN(p, :, 1));
        dN_deta = squeeze(dN(p, :, 2));
        
        if size(dN_dxi, 1) > size(dN_dxi, 2)
            dN_dxi = dN_dxi';
        end
        if size(dN_deta, 1) > size(dN_deta, 2)
            dN_deta = dN_deta';
        end
        
        J(1, 1, p) = dN_dxi * nodeCoords(:, 1);
        J(1, 2, p) = dN_deta * nodeCoords(:, 1);
        J(2, 1, p) = dN_dxi * nodeCoords(:, 2);
        J(2, 2, p) = dN_deta * nodeCoords(:, 2);
        
        detJ(p) = det(J(:, :, p));
    end
end

if n_points == 1
    J = J(:, :, 1);
    detJ = detJ(1);
end

end

