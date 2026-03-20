function [K, M] = assembleGlobalMatrices(nodes, elements, N)
% assembleGlobalMatrices - Assemble global stiffness and mass matrices
%
% Input:
%   nodes    - Node data structure
%   elements - Element data structure
%   N        - Number structure (contains N.n, N.el, etc.)
%
% Output:
%   K - Global stiffness matrix [n_dof × n_dof] (sparse)
%   M - Global mass matrix [n_dof × n_dof] (sparse)

n_dof = N.n;

K = sparse(n_dof, n_dof);
M = sparse(n_dof, n_dof);

for el = 1:N.el
    el_nodes = elements.all(el, :);
    el_nodes = el_nodes(el_nodes ~= 0);
    
    node_indices = nodes.nind(el_nodes, 2);
    node_coords = nodes.all(node_indices, 1:2);
    
    element_type = elements.type(el);
    
    [K_elem, M_elem] = computeElementMatrices(node_coords, element_type);
    
    for i = 1:length(el_nodes)
        for j = 1:length(el_nodes)
            I = nodes.nind(el_nodes(i), 2);
            J = nodes.nind(el_nodes(j), 2);
            
            K(I, J) = K(I, J) + K_elem(i, j);
            M(I, J) = M(I, J) + M_elem(i, j);
        end
    end
end

end

