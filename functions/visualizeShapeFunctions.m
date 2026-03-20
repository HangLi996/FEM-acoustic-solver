function visualizeShapeFunctions(elementType, nPoints)
% visualizeShapeFunctions - Visualize Lagrangian shape functions
%
% Input:
%   elementType - Element type: 1 (Line2), 8 (Line3), 3 (Quadrangle4), 10 (Quadrangle9)
%   nPoints     - Optional, number of grid points for plotting (default: 50)
%
% Output:
%   None, displays figure window directly
%
% Function:
%   Visualize shape functions on reference coordinate system
%   - 1D elements: plot shape function curves
%   - 2D elements: plot shape function 3D surfaces

if nargin < 2
    nPoints = 50;  % Default grid points
end

switch elementType
    case 1  % Line2 - first-order line element (2 nodes)
        visualizeLine2(nPoints);
        
    case 8  % Line3 - second-order line element (3 nodes)
        visualizeLine3(nPoints);
        
    case 3  % Quadrangle4 - first-order quadrilateral (4 nodes)
        visualizeQuad4(nPoints);
        
    case 10 % Quadrangle9 - second-order quadrilateral (9 nodes)
        visualizeQuad9(nPoints);
        
    otherwise
        error('Unsupported element type: %d', elementType);
end

end

%% Line2 visualization
function visualizeLine2(nPoints)
xi = linspace(-1, 1, nPoints);
[N, ~] = shapeFunctions(xi, [], 1);

figure('Name', 'Line2 Shape Functions', 'Position', [100, 100, 1200, 400]);

% Plot shape functions
subplot(1, 2, 1);
plot(xi, N(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(xi, N(:, 2), 'r-', 'LineWidth', 2);
plot([-1, -1], [0, 1], 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot([1, 1], [0, 1], 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('ξ');
ylabel('N_i(ξ)');
title('Line2 - Shape Functions');
legend('N_0', 'N_1', 'Location', 'best');
grid on;

% Plot sum of shape functions (verify partition of unity)
subplot(1, 2, 2);
plot(xi, sum(N, 2), 'g-', 'LineWidth', 2);
xlabel('ξ');
ylabel('Σ N_i(ξ)');
title('Partition of Unity Verification (should be 1)');
ylim([0.99, 1.01]);
grid on;
end

%% Line3 visualization
function visualizeLine3(nPoints)
xi = linspace(-1, 1, nPoints);
[N, ~] = shapeFunctions(xi, [], 8);

figure('Name', 'Line3 Shape Functions', 'Position', [100, 100, 1200, 400]);

% Plot shape functions
subplot(1, 2, 1);
plot(xi, N(:, 1), 'b-', 'LineWidth', 2); hold on;
plot(xi, N(:, 2), 'g-', 'LineWidth', 2);
plot(xi, N(:, 3), 'r-', 'LineWidth', 2);
plot([-1, -1], [0, 1], 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot([0, 0], [0, 1], 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot([1, 1], [0, 1], 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('ξ');
ylabel('N_i(ξ)');
title('Line3 - Shape Functions');
legend('N_0', 'N_2', 'N_1', 'Location', 'best');
grid on;

% Plot sum of shape functions
subplot(1, 2, 2);
plot(xi, sum(N, 2), 'g-', 'LineWidth', 2);
xlabel('ξ');
ylabel('Σ N_i(ξ)');
title('Partition of Unity Verification (should be 1)');
ylim([0.99, 1.01]);
grid on;
end

%% Quadrangle4 visualization
function visualizeQuad4(nPoints)
xi = linspace(-1, 1, nPoints);
eta = linspace(-1, 1, nPoints);
[XI, ETA] = meshgrid(xi, eta);

% Compute shape function values
N = zeros(nPoints, nPoints, 4);
for i = 1:nPoints
    for j = 1:nPoints
        [N(i, j, :), ~] = shapeFunctions(XI(i, j), ETA(i, j), 3);
    end
end

figure('Name', 'Quadrangle4 Shape Functions', 'Position', [50, 50, 1400, 700]);

% Node coordinates
node_coords = [-1 -1; 1 -1; 1 1; -1 1];

% Plot 4 shape functions
for i = 1:4
    subplot(2, 2, i);
    surf(XI, ETA, N(:, :, i));
    hold on;
    % Mark node position
    plot3(node_coords(i, 1), node_coords(i, 2), 1, 'ro', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r');
    % Mark other nodes (value = 0)
    other_nodes = setdiff(1:4, i);
    plot3(node_coords(other_nodes, 1), node_coords(other_nodes, 2), ...
        zeros(size(other_nodes)), 'ko', 'MarkerSize', 8);
    
    xlabel('ξ');
    ylabel('η');
    zlabel(sprintf('N_%d(ξ,η)', i-1));
    title(sprintf('Quadrangle4 - N_%d', i-1));
    shading interp;
    colormap(jet);
    colorbar;
    view(45, 30);
end
end

%% Quadrangle9 visualization
function visualizeQuad9(nPoints)
xi = linspace(-1, 1, nPoints);
eta = linspace(-1, 1, nPoints);
[XI, ETA] = meshgrid(xi, eta);

% Compute shape function values
N = zeros(nPoints, nPoints, 9);
for i = 1:nPoints
    for j = 1:nPoints
        [N(i, j, :), ~] = shapeFunctions(XI(i, j), ETA(i, j), 10);
    end
end

figure('Name', 'Quadrangle9 Shape Functions', 'Position', [20, 20, 1600, 1200]);

% Node coordinates (GMSH convention)
% Corner nodes: 0(-1,-1), 1(1,-1), 2(1,1), 3(-1,1)
% Mid-side nodes: 4(0,-1), 5(1,0), 6(0,1), 7(-1,0)
% Center: 8(0,0)
node_coords = [-1 -1; 1 -1; 1 1; -1 1; ...  % Corner nodes 0-3
               0 -1; 1 0; 0 1; -1 0; ...   % Mid-side nodes 4-7
               0 0];                        % Center 8

% Plot 9 shape functions (3x3 grid)
for i = 1:9
    subplot(3, 3, i);
    surf(XI, ETA, N(:, :, i));
    hold on;
    
    % Mark current node position (value = 1)
    plot3(node_coords(i, 1), node_coords(i, 2), 1, 'ro', ...
        'MarkerSize', 12, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    
    % Mark other node positions (value = 0)
    other_nodes = setdiff(1:9, i);
    plot3(node_coords(other_nodes, 1), node_coords(other_nodes, 2), ...
        zeros(size(other_nodes)), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
    
    xlabel('ξ');
    ylabel('η');
    zlabel(sprintf('N_%d(ξ,η)', i-1));
    title(sprintf('N_%d', i-1), 'FontSize', 12, 'FontWeight', 'bold');
    shading interp;
    colormap(jet);
    colorbar;
    view(45, 30);
    zlim([0, 1.1]);
end

% Add overall title
sgtitle('Second-order Lagrangian Shape Functions N_i (i=0,...,8) for Quadrangle9', ...
    'FontSize', 14, 'FontWeight', 'bold');
end
