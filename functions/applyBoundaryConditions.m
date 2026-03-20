function [K_mod, M_mod, f_mod] = applyBoundaryConditions(K, M, f, boundary_nodes, boundary_type, nodes)
% applyBoundaryConditions - Apply boundary conditions to system matrices
%
% Input:
%   K              - Global stiffness matrix [n_dof × n_dof]
%   M              - Global mass matrix [n_dof × n_dof]
%   f              - Load vector [n_dof × 1] (optional)
%   boundary_nodes - Boundary node indices (global node numbers)
%   boundary_type  - Boundary type: 'hard' (Dirichlet, p=0) or 'soft' (Neumann, natural)
%   nodes          - Node data structure (for node index mapping)
%
% Output:
%   K_mod - Modified stiffness matrix
%   M_mod - Modified mass matrix
%   f_mod - Modified load vector

if nargin < 6
    error('Need nodes structure for index mapping');
end

if nargin < 5
    boundary_type = 'hard';
end

% Normalize boundary_nodes to numeric vector
if iscell(boundary_nodes)
    boundary_nodes = unique(cell2mat(boundary_nodes(:)));
end
boundary_nodes = boundary_nodes(:);

if nargin < 4 || isempty(boundary_nodes)
    K_mod = K;
    M_mod = M;
    if nargin >= 3 && ~isempty(f)
        f_mod = f;
    else
        f_mod = [];
    end
    return;
end

K_mod = K;
M_mod = M;

if nargin >= 3 && ~isempty(f)
    f_mod = f;
else
    f_mod = sparse(size(K, 1), 1);
end

% Convert boundary node numbers to DOF indices
boundary_dof = nodes.nind(boundary_nodes, 2);
boundary_dof = boundary_dof(boundary_dof > 0);

if strcmpi(boundary_type, 'hard') || strcmpi(boundary_type, 'dirichlet')
    % Dirichlet boundary condition: p = 0
    for idx = boundary_dof(:)'
        % Modify stiffness matrix
        K_mod(idx, :) = 0;
        K_mod(:, idx) = 0;
        K_mod(idx, idx) = 1;
        
        % Modify mass matrix
        M_mod(idx, :) = 0;
        M_mod(:, idx) = 0;
        M_mod(idx, idx) = 0;
        
        % Modify load vector
        f_mod(idx) = 0;
    end
elseif strcmpi(boundary_type, 'soft') || strcmpi(boundary_type, 'neumann')
    % Neumann boundary condition: natural boundary condition
    % No modification needed (automatically satisfied in weak form)
    % Just return original matrices
else
    warning('Unknown boundary type: %s. Using natural boundary condition.', boundary_type);
end

end

