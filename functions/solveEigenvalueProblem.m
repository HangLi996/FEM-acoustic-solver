function [eigenvalues, eigenvectors, eigenfrequencies] = solveEigenvalueProblem(K, M, air, n_modes)
% solveEigenvalueProblem - Solve acoustic eigenvalue problem
%
% Input:
%   K        - Stiffness matrix [n_dof × n_dof] (sparse or full)
%   M        - Mass matrix [n_dof × n_dof] (sparse or full)
%   air      - Air parameters structure (contains speed of sound c)
%   n_modes  - Optional, number of modes to compute (default: 10)
%
% Output:
%   eigenvalues      - Eigenvalues [n_modes × 1] (k² values)
%   eigenvectors     - Eigenvectors [n_dof × n_modes] (pressure modes)
%   eigenfrequencies  - Eigenfrequencies [n_modes × 1] (Hz)

if nargin < 4
    n_modes = 10;
end

if nargin < 3 || ~isfield(air, 'c')
    error('Need air structure with speed of sound c');
end

if size(K, 1) ~= size(K, 2) || size(M, 1) ~= size(M, 2)
    error('K and M must be square matrices');
end

if size(K, 1) ~= size(M, 1)
    error('K and M must have same dimensions');
end

n_dof = size(K, 1);

n_modes = min(n_modes, n_dof);

try
    [eigenvectors, eigenvals] = eigs(K, M, n_modes, 'smallestabs');
catch ME
    warning('eigs failed, trying eig (may be slower)');
    [eigenvectors, eigenvals] = eig(full(K), full(M));
    [~, idx] = sort(diag(eigenvals));
    idx = idx(1:min(n_modes, length(idx)));
    eigenvals = eigenvals(idx, idx);
    eigenvectors = eigenvectors(:, idx);
end

eigenvalues = diag(eigenvals);

eigenvalues = real(eigenvalues);
eigenvalues(eigenvalues < 0) = NaN;

k = sqrt(eigenvalues);

eigenfrequencies = k * air.c / (2 * pi);

[eigenfrequencies, idx] = sort(eigenfrequencies);
eigenvalues = eigenvalues(idx);
eigenvectors = eigenvectors(:, idx);

valid_idx = ~isnan(eigenfrequencies);
eigenvalues = eigenvalues(valid_idx);
eigenvectors = eigenvectors(:, valid_idx);
eigenfrequencies = eigenfrequencies(valid_idx);

end

