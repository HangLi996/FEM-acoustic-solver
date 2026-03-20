function p = solveLinearSystem(A, f, method)
% solveLinearSystem - Solve linear system A·p = f
%
% Input:
%   A      - System matrix [n_dof × n_dof]
%   f      - Right-hand side vector [n_dof × 1]
%   method - Solution method: 'direct' (default) or 'iterative'
%
% Output:
%   p - Solution vector [n_dof × 1]

if nargin < 3
    method = 'direct';
end

if size(A, 1) ~= size(A, 2)
    error('System matrix A must be square');
end

if size(A, 1) ~= size(f, 1)
    error('Dimensions of A and f must match');
end

if strcmpi(method, 'direct')
    % Direct method using MATLAB backslash (automatic solver selection)
    p = A \ f;
    
elseif strcmpi(method, 'iterative')
    % Iterative method (GMRES for general matrices)
    tol = 1e-6;
    maxit = min(1000, size(A, 1));
    
    try
        [p, flag, relres, iter] = gmres(A, f, [], tol, maxit);
        if flag ~= 0
            warning('GMRES did not converge: flag=%d, relres=%.2e', flag, relres);
        end
    catch ME
        warning('GMRES failed, using direct method: %s', ME.message);
        p = A \ f;
    end
    
else
    error('Unknown solution method: %s', method);
end

end

