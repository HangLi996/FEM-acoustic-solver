function convergenceStudy(mesh_files, element_types, air, n_modes)
% convergenceStudy - Study eigenvalue convergence with mesh refinement
%
% Input:
%   mesh_files    - List of mesh files (coarse to fine) {'mesh1.msh', 'mesh2.msh', ...}
%   element_types - List of element types [type1, type2, ...]
%   air           - Air parameters structure (contains speed of sound c)
%   n_modes       - Number of modes to study (default: 5)
%
% Output:
%   None, displays convergence plots

if nargin < 4
    n_modes = 5;
end

if nargin < 3 || ~isfield(air, 'c')
    error('Need air structure with speed of sound c');
end

n_meshes = length(mesh_files);
n_types = length(element_types);

eigenvalues_history = cell(n_types, n_meshes);
element_sizes = zeros(n_meshes, 1);
h_labels = cell(n_meshes, 1);

for type_idx = 1:n_types
    element_type = element_types(type_idx);
    
    fprintf('=== Element type %d ===\n', element_type);
    
    for mesh_idx = 1:n_meshes
        fprintf('Processing mesh %d/%d: %s\n', mesh_idx, n_meshes, mesh_files{mesh_idx});
        
        file.mesh = mesh_files{mesh_idx};
        if ~isfield(file, 'path')
            file.path = './';
        end
        conf.mesher = 'GMSH';
        [nodes, elements, dim, N] = msh2mat(file);
        [B, N] = getBoundary(file, conf, N, nodes, elements);
        
        h_sum = 0;
        h_count = 0;
        for el = 1:N.el
            el_nodes = elements.all(el, :);
            el_nodes = el_nodes(el_nodes ~= 0);
            node_coords = nodes.all(nodes.nind(el_nodes, 2), 1:2);
            
            if length(el_nodes) >= 2
                h = norm(node_coords(end, :) - node_coords(1, :));
                h_sum = h_sum + h;
                h_count = h_count + 1;
            end
        end
        element_sizes(mesh_idx) = h_sum / h_count;
        h_labels{mesh_idx} = sprintf('h=%.3f', element_sizes(mesh_idx));
        
        [K, M] = assembleGlobalMatrices(nodes, elements, N);
        
        [eigenvalues, ~, eigenfrequencies] = solveEigenvalueProblem(K, M, air, n_modes);
        
        eigenvalues_history{type_idx, mesh_idx} = eigenfrequencies;
    end
end

figure;
set(gcf, 'Position', [100, 100, 1400, 600]);

for type_idx = 1:n_types
    subplot(1, n_types, type_idx);
    hold on;
    
    freq_data = zeros(n_modes, n_meshes);
    for mesh_idx = 1:n_meshes
        freqs = eigenvalues_history{type_idx, mesh_idx};
        n_freqs = min(length(freqs), n_modes);
        freq_data(1:n_freqs, mesh_idx) = freqs(1:n_freqs);
    end
    
    colors = lines(n_modes);
    for mode = 1:n_modes
        plot(element_sizes, freq_data(mode, :), 'o-', ...
            'Color', colors(mode, :), 'LineWidth', 2, ...
            'MarkerSize', 8, 'DisplayName', sprintf('Mode %d', mode));
    end
    
    xlabel('Element size h (m)', 'FontSize', 12);
    ylabel('Eigenfrequency f (Hz)', 'FontSize', 12);
    title(sprintf('Convergence for element type %d', element_types(type_idx)), 'FontSize', 14);
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 11);
end

sgtitle('Eigenvalue convergence study', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('\n=== Convergence study summary ===\n');
for type_idx = 1:n_types
    fprintf('\nElement type %d:\n', element_types(type_idx));
    fprintf('%-15s', 'Mesh');
    for mode = 1:min(5, n_modes)
        fprintf('  Mode%d(Hz)', mode);
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for mesh_idx = 1:n_meshes
        freqs = eigenvalues_history{type_idx, mesh_idx};
        fprintf('%-15s', h_labels{mesh_idx});
        for mode = 1:min(5, n_modes)
            if mode <= length(freqs)
                fprintf('  %8.2f', freqs(mode));
            else
                fprintf('  %8s', 'N/A');
            end
        end
        fprintf('\n');
    end
end

end

