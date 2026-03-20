function [K_elem, M_elem] = computeElementMatrices(nodeCoords, elementType)
% computeElementMatrices - Compute element stiffness and mass matrices
%
% Input:
%   nodeCoords  - Element node coordinates [n_nodes × 2]
%   elementType - Element type: 1 (Line2), 8 (Line3), 3 (Quadrangle4), 10 (Quadrangle9)
%
% Output:
%   K_elem - Element stiffness matrix [n_nodes × n_nodes]
%   M_elem - Element mass matrix [n_nodes × n_nodes]

is1D = (elementType == 1 || elementType == 8);
is2D = (elementType == 3 || elementType == 10);

if ~is1D && ~is2D
    error('Unsupported element type: %d', elementType);
end

n_nodes = size(nodeCoords, 1);

K_elem = zeros(n_nodes, n_nodes);
M_elem = zeros(n_nodes, n_nodes);

% Gauss quadrature points and weights
if is1D
    gauss_points = [-1/sqrt(3), 1/sqrt(3)];
    gauss_weights = [1, 1];
    
    for i = 1:length(gauss_points)
        xi = gauss_points(i);
        w_i = gauss_weights(i);
        
        [N, dN] = shapeFunctions(xi, [], elementType);
        dN_dxi = squeeze(dN(1, :, 1));
        if size(dN_dxi, 1) == 1
            dN_dxi = dN_dxi';
        end
        
        [J, detJ] = computeJacobian(xi, [], nodeCoords, elementType);
        dN_dx = dN_dxi / J;
        weight = w_i * abs(detJ);
        
        for a = 1:n_nodes
            for b = 1:n_nodes
                K_elem(a, b) = K_elem(a, b) + dN_dx(a) * dN_dx(b) * weight;
                M_elem(a, b) = M_elem(a, b) + N(1, a) * N(1, b) * weight;
            end
        end
    end
    
else
    gauss_points = [-1/sqrt(3), 1/sqrt(3)];
    gauss_weights = [1, 1];
    
    for i = 1:length(gauss_points)
        for j = 1:length(gauss_points)
            xi = gauss_points(i);
            eta = gauss_points(j);
            w_i = gauss_weights(i);
            w_j = gauss_weights(j);
            
            [N, dN] = shapeFunctions(xi, eta, elementType);
            dN_dxi = squeeze(dN(1, :, 1));
            dN_deta = squeeze(dN(1, :, 2));
            
            if size(dN_dxi, 1) > size(dN_dxi, 2)
                dN_dxi = dN_dxi';
            end
            if size(dN_deta, 1) > size(dN_deta, 2)
                dN_deta = dN_deta';
            end
            
            [J, detJ] = computeJacobian(xi, eta, nodeCoords, elementType);
            J_inv = inv(J);
            
            dN_dx = dN_dxi * J_inv(1, 1) + dN_deta * J_inv(1, 2);
            dN_dy = dN_dxi * J_inv(2, 1) + dN_deta * J_inv(2, 2);
            
            weight = w_i * w_j * abs(detJ);
            
            for a = 1:n_nodes
                for b = 1:n_nodes
                    K_elem(a, b) = K_elem(a, b) + ...
                        (dN_dx(a) * dN_dx(b) + dN_dy(a) * dN_dy(b)) * weight;
                    M_elem(a, b) = M_elem(a, b) + ...
                        N(1, a) * N(1, b) * weight;
                end
            end
        end
    end
end

end

