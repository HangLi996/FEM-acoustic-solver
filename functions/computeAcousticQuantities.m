function [SPL, I, W] = computeAcousticQuantities(p, nodes, elements, B, omega, air)
% computeAcousticQuantities - Compute acoustic quantities from pressure field
%
% Input:
%   p      - Pressure field [n_dof × 1] or [n_dof × n_freq]
%   nodes  - Node data structure
%   elements - Element data structure
%   B      - Boundary structure (optional)
%   omega  - Angular frequency [rad/s] or vector [n_freq × 1]
%   air    - Air parameters structure (contains rho, c)
%
% Output:
%   SPL - Sound pressure level [n_dof × 1] or [n_dof × n_freq] (dB)
%   I   - Sound intensity (optional, requires velocity field)
%   W   - Sound power (optional, requires boundary integration)

if nargin < 6 || ~isfield(air, 'rho') || ~isfield(air, 'c')
    error('Need air structure with rho and c');
end

% Reference pressure
p_ref = 2e-5;  % [Pa]

% Sound pressure level
SPL = 20 * log10(abs(p) / p_ref);

% Sound intensity (simplified, requires velocity field for full calculation)
I = [];
if nargout >= 2
    % Simplified: I = |p|² / (2*rho*c) for plane waves
    I = abs(p).^2 / (2 * air.rho * air.c);
end

% Sound power (requires boundary integration)
W = [];
if nargout >= 3 && nargin >= 4 && ~isempty(B)
    % Simplified calculation on boundary
    % Full implementation would require velocity field and boundary integration
    warning('Sound power calculation requires velocity field - returning empty');
end

end

