function [N, dN] = shapeFunctions(xi, eta, elementType)
% shapeFunctions - Compute Lagrangian shape functions and their gradients
%
% Input:
%   xi          - Natural coordinate xi (scalar or vector)
%   eta         - Natural coordinate eta (scalar or vector, for 2D elements only)
%   elementType - Element type: 1 (Line2), 8 (Line3), 3 (Quadrangle4), 10 (Quadrangle9)
%
% Output:
%   N  - Shape function values [n_points × n_nodes]
%   dN - Shape function gradients [n_points × n_nodes × n_dim]
%        For 1D: dN(:, :, 1) = dN/dxi
%        For 2D: dN(:, :, 1) = dN/dxi, dN(:, :, 2) = dN/deta

% Convert to column vectors
xi = xi(:);
if nargin >= 2 && ~isempty(eta)
    eta = eta(:);
end

switch elementType
    case 1  % Line2 - first-order line element (2 nodes)
        [N, dN] = shapeLine2(xi);
        
    case 8  % Line3 - second-order line element (3 nodes)
        [N, dN] = shapeLine3(xi);
        
    case 3  % Quadrangle4 - first-order quadrilateral (4 nodes)
        if nargin < 2 || isempty(eta)
            error('Quadrangle4 requires eta parameter');
        end
        [N, dN] = shapeQuad4(xi, eta);
        
    case 10 % Quadrangle9 - second-order quadrilateral (9 nodes)
        if nargin < 2 || isempty(eta)
            error('Quadrangle9 requires eta parameter');
        end
        [N, dN] = shapeQuad9(xi, eta);
        
    otherwise
        error('Unsupported element type: %d', elementType);
end

end

%% Line2 - first-order line element (2 nodes)
function [N, dN] = shapeLine2(xi)
% Node numbering: 0 (xi=-1), 1 (xi=1)
% Shape functions:
%   N0 = (1 - xi)/2
%   N1 = (1 + xi)/2

n_points = length(xi);
N = zeros(n_points, 2);
dN = zeros(n_points, 2, 1);

% Node 0 (xi = -1)
N(:, 1) = (1 - xi) / 2;
dN(:, 1, 1) = -1/2;

% Node 1 (xi = 1)
N(:, 2) = (1 + xi) / 2;
dN(:, 2, 1) = 1/2;

end

%% Line3 - second-order line element (3 nodes)
function [N, dN] = shapeLine3(xi)
% Node numbering: 0 (xi=-1), 2 (xi=0), 1 (xi=1)
% Shape functions:
%   N0 = xi(xi - 1)/2
%   N2 = (1 - xi^2)
%   N1 = xi(xi + 1)/2

n_points = length(xi);
N = zeros(n_points, 3);
dN = zeros(n_points, 3, 1);

% Node 0 (xi = -1)
N(:, 1) = xi .* (xi - 1) / 2;
dN(:, 1, 1) = (2*xi - 1) / 2;

% Node 2 (xi = 0)
N(:, 2) = 1 - xi.^2;
dN(:, 2, 1) = -2*xi;

% Node 1 (xi = 1)
N(:, 3) = xi .* (xi + 1) / 2;
dN(:, 3, 1) = (2*xi + 1) / 2;

end

%% Quadrangle4 - first-order quadrilateral (4 nodes)
function [N, dN] = shapeQuad4(xi, eta)
% Node numbering:
%   0: (-1, -1)  1: (1, -1)
%   3: (-1, 1)   2: (1, 1)
%
% Shape functions (bilinear):
%   N0 = (1 - xi)(1 - eta)/4
%   N1 = (1 + xi)(1 - eta)/4
%   N2 = (1 + xi)(1 + eta)/4
%   N3 = (1 - xi)(1 + eta)/4

n_points = length(xi);
N = zeros(n_points, 4);
dN = zeros(n_points, 4, 2);

% Node 0: (-1, -1)
N(:, 1) = (1 - xi) .* (1 - eta) / 4;
dN(:, 1, 1) = -(1 - eta) / 4;  % dN/dxi
dN(:, 1, 2) = -(1 - xi) / 4;    % dN/deta

% Node 1: (1, -1)
N(:, 2) = (1 + xi) .* (1 - eta) / 4;
dN(:, 2, 1) = (1 - eta) / 4;    % dN/dxi
dN(:, 2, 2) = -(1 + xi) / 4;    % dN/deta

% Node 2: (1, 1)
N(:, 3) = (1 + xi) .* (1 + eta) / 4;
dN(:, 3, 1) = (1 + eta) / 4;    % dN/dxi
dN(:, 3, 2) = (1 + xi) / 4;     % dN/deta

% Node 3: (-1, 1)
N(:, 4) = (1 - xi) .* (1 + eta) / 4;
dN(:, 4, 1) = -(1 + eta) / 4;   % dN/dxi
dN(:, 4, 2) = (1 - xi) / 4;     % dN/deta

end

%% Quadrangle9 - second-order quadrilateral (9 nodes)
function [N, dN] = shapeQuad9(xi, eta)
% Node number:
%   0: (-1, -1)  4: (0, -1)   1: (1, -1)
%   7: (-1, 0)   8: (0, 0)    5: (1, 0)
%   3: (-1, 1)   6: (0, 1)    2: (1, 1)
%
% Shape functions (biquadratic):
%   Corner nodes (i=0,1,2,3):
%     Ni = (1/4) * xi_i*xi * (1 + xi_i*xi) * eta_i*eta * (1 + eta_i*eta)
%   Mid-side nodes (i=4,5,6,7):
%     N4 = (1 - xi^2) * (-eta(1-eta)/2)  [bottom edge]
%     N5 = (xi(1+xi)/2) * (1 - eta^2)   [right edge]
%     N6 = (1 - xi^2) * (eta(1+eta)/2)   [top edge]
%     N7 = (xi(xi-1)/2) * (1 - eta^2)  [left edge]
%   Center node (i=8):
%     N8 = (1 - xi^2) * (1 - eta^2)

n_points = length(xi);
N = zeros(n_points, 9);
dN = zeros(n_points, 9, 2);

% Corner node coordinates
xi_corner = [-1, 1, 1, -1];
eta_corner = [-1, -1, 1, 1];

% Corner nodes (0, 1, 2, 3)
for i = 1:4
    xi_i = xi_corner(i);
    eta_i = eta_corner(i);
    N(:, i) = (1/4) * (xi_i.*xi) .* (1 + xi_i.*xi) .* (eta_i.*eta) .* (1 + eta_i.*eta);
    
    dN(:, i, 1) = (1/4) * xi_i .* (1 + 2*xi_i.*xi) .* (eta_i.*eta) .* (1 + eta_i.*eta);
    dN(:, i, 2) = (1/4) * (xi_i.*xi) .* (1 + xi_i.*xi) .* eta_i .* (1 + 2*eta_i.*eta);
end

% Mid-side nodes
% Node 4: bottom edge midpoint (0, -1)
N(:, 5) = (1 - xi.^2) .* (-eta .* (1 - eta) / 2);
dN(:, 5, 1) = -2*xi .* (-eta .* (1 - eta) / 2);
dN(:, 5, 2) = (1 - xi.^2) .* (-(1 - 2*eta) / 2);

% Node 5: right edge midpoint (1, 0)
N(:, 6) = (xi .* (1 + xi) / 2) .* (1 - eta.^2);
dN(:, 6, 1) = ((1 + 2*xi) / 2) .* (1 - eta.^2);
dN(:, 6, 2) = (xi .* (1 + xi) / 2) .* (-2*eta);

% Node 6: top edge midpoint (0, 1)
N(:, 7) = (1/2) * (1 - xi.^2) .* eta .* (1 + eta);
dN(:, 7, 1) = -xi .* eta .* (1 + eta);
dN(:, 7, 2) = (1/2) * (1 - xi.^2) .* (1 + 2*eta);

% Node 7: left edge midpoint (-1, 0)
N(:, 8) = (xi .* (xi - 1) / 2) .* (1 - eta.^2);
dN(:, 8, 1) = ((2*xi - 1) / 2) .* (1 - eta.^2);
dN(:, 8, 2) = (xi .* (xi - 1) / 2) .* (-2*eta);

% Center node 8: (0, 0)
N(:, 9) = (1 - xi.^2) .* (1 - eta.^2);
dN(:, 9, 1) = -2*xi .* (1 - eta.^2);
dN(:, 9, 2) = -2*eta .* (1 - xi.^2);

end
