% run_frequency_response.m — 声学有限元：频响与谐波声场示例
%
% 演示：施加边界条件、扫频频响、在选定频率下求解声压场并可视化。
% 对应离散 Helmholtz 系统 (K - k²M)p = f，详见 functions/ 中文档说明。
%
% 用法：在仓库根目录执行 run_frequency_response
%
% 依赖：functions/；默认网格见下方 file.mesh

clear all; close all; clc;
addpath('functions');
addpath('mesh');
config;

%% 用户可改：网格文件
file.mesh = 'mesh/duct/duct4l.msh';

fprintf('========================================\n');
fprintf('  FEM 声学求解器 — 频响与声场示例\n');
fprintf('========================================\n\n');

%% 加载网格与边界
fprintf('=== 加载网格 ===\n');
[nodes, elements, dim, N] = msh2mat(file);
[B, N] = getBoundary(file, conf, N, nodes, elements);

fprintf('网格信息：\n');
fprintf('  节点数：%d\n', N.n);
fprintf('  单元数：%d\n', N.el);
fprintf('  空间维数：%d\n', dim);
fprintf('  边界节点数：%d\n', length(B.n));
fprintf('\n');

%% 组装全局矩阵
fprintf('=== 组装全局矩阵 ===\n');
fprintf('正在计算单元矩阵并组装…\n');
tic;
[K, M] = assembleGlobalMatrices(nodes, elements, N);
time_assembly = toc;
fprintf('完成（耗时 %.4f 秒）\n', time_assembly);
fprintf('  全局刚度矩阵 K：[%d × %d]\n', size(K, 1), size(K, 2));
fprintf('  全局质量矩阵 M：[%d × %d]\n', size(M, 1), size(M, 2));
fprintf('\n');

%% 施加边界条件（硬边界：Dirichlet p=0）
fprintf('========================================\n');
fprintf('步骤 1：边界条件\n');
fprintf('========================================\n');

fprintf('施加硬边界（Dirichlet，p = 0）…\n');
f_init = sparse(N.n, 1);
[K, M, f_init] = applyBoundaryConditions(K, M, f_init, B.n, 'hard', nodes);
fprintf('完成。\n');
if iscell(B.n)
    boundary_count = sum(cellfun(@length, B.n));
    fprintf('  边界节点：%d（来自 %d 条边界）\n', boundary_count, length(B.n));
else
    fprintf('  边界节点：%d\n', length(B.n));
end
fprintf('  边界类型：硬边界（p = 0）\n');
fprintf('\n');

%% 频响分析
fprintf('========================================\n');
fprintf('步骤 2：频响分析\n');
fprintf('========================================\n');

fprintf('添加点声源…\n');
desired_source_position = [1.7, 0];  % 目标位置（约 3.4 m 管道中部）

if iscell(B.n)
    boundary_node_list = unique(cell2mat(B.n(:)));
else
    boundary_node_list = B.n(:);
end
all_nodes = 1:size(nodes.all, 1);
interior_nodes = setdiff(all_nodes, boundary_node_list);

node_coords = nodes.all(:, 1:2);
interior_coords = node_coords(interior_nodes, :);
distances = sqrt(sum((interior_coords - repmat(desired_source_position(1:2), size(interior_coords, 1), 1)).^2, 2));
[~, closest_interior_idx] = min(distances);
source_node = interior_nodes(closest_interior_idx);
source_position = node_coords(source_node, 1:2);

source_strength = 1.0;
fprintf('  目标位置：[%.2f, %.2f] m\n', desired_source_position(1), desired_source_position(2));
fprintf('  实际声源（最近内部节点）：[%.2f, %.2f] m（节点 %d）\n', ...
    source_position(1), source_position(2), source_node);
fprintf('  声源强度 Q：%.2f\n', source_strength);

freq_range = [10, 500, 200];  % [f_min, f_max, 采样点数]
fprintf('正在计算频响…\n');
fprintf('  频率范围：%.0f – %.0f Hz\n', freq_range(1), freq_range(2));
fprintf('  采样点数：%d\n', freq_range(3));

source_info.nodes = nodes;
source_info.source_position = source_position;
source_info.source_strength = source_strength;

[p_freq, frequencies] = computeFrequencyResponse(K, M, source_info, freq_range, air, nodes, B.n, 'hard');
fprintf('完成。\n\n');

source_node = findClosestNode(nodes, source_position);
source_dof = nodes.nind(source_node, 2);
if iscell(B.n)
    boundary_node_list = unique(cell2mat(B.n(:)));
else
    boundary_node_list = B.n(:);
end
is_source_on_boundary = ismember(source_node, boundary_node_list);

if is_source_on_boundary
    warning('声源落在边界节点上，响应可能被抑制。');
    fprintf('  声源节点 %d 为边界节点。\n', source_node);
else
    fprintf('  声源节点 %d 不在边界上。\n', source_node);
end

obs_position = [0.85, 0];  % 观测点（约 1/4 管长）
obs_node = findClosestNode(nodes, obs_position);
obs_dof = nodes.nind(obs_node, 2);

is_obs_on_boundary = ismember(obs_node, boundary_node_list);
if is_obs_on_boundary
    fprintf('  注意：观测点节点 %d 在边界上，改用内部最近节点。\n', obs_node);
    interior_coords = node_coords(interior_nodes, :);
    obs_distances = sqrt(sum((interior_coords - repmat(obs_position(1:2), size(interior_coords, 1), 1)).^2, 2));
    [~, closest_interior_obs_idx] = min(obs_distances);
    obs_node = interior_nodes(closest_interior_obs_idx);
    obs_dof = nodes.nind(obs_node, 2);
    fprintf('  改用内部观测点：节点 %d，[%.2f, %.2f] m\n', ...
        obs_node, nodes.all(obs_node, 1), nodes.all(obs_node, 2));
end

p_freq_full = full(p_freq);
fprintf('  p_freq 尺寸：[%d × %d]\n', size(p_freq, 1), size(p_freq, 2));
fprintf('  幅值统计：min=%.6e, max=%.6e, mean(|·|)=%.6e\n', ...
    min(p_freq_full(:)), max(p_freq_full(:)), mean(abs(p_freq_full(:))));

if obs_dof > 0 && obs_dof <= size(p_freq, 1)
    response_raw = abs(p_freq(obs_dof, :));
    fprintf('  观测点响应范围：[%.6e, %.6e] Pa\n', min(response_raw), max(response_raw));
    fprintf('  观测点：节点 %d（DOF %d），[%.2f, %.2f] m\n', ...
        obs_node, obs_dof, nodes.all(obs_node, 1), nodes.all(obs_node, 2));
else
    warning('观测 DOF %d 无效（最大 %d），改用全体平均响应。', obs_dof, size(p_freq, 1));
    response_raw = mean(abs(p_freq), 1);
end

fprintf('  声源 DOF：%d\n', source_dof);

omega_vec = 2 * pi * frequencies;
response = response_raw ./ omega_vec;
fprintf('  归一化响应范围：[%.6e, %.6e] Pa·s/rad\n', ...
    min(response), max(response));

if iscell(B.n)
    boundary_node_list = unique(cell2mat(B.n(:)));
else
    boundary_node_list = B.n(:);
end
all_nodes = 1:size(nodes.all, 1);
interior_nodes = setdiff(all_nodes, boundary_node_list);
interior_dof = nodes.nind(interior_nodes, 2);
interior_dof = interior_dof(interior_dof > 0);

if ~isempty(interior_dof)
    valid_dof = interior_dof(interior_dof > 0 & interior_dof <= size(p_freq, 1));
    if ~isempty(valid_dof)
        response_avg = mean(abs(p_freq(valid_dof, :)), 1);
        fprintf('  使用 %d 个内部节点的平均响应。\n', length(valid_dof));
        fprintf('  平均响应范围：[%.6e, %.6e] Pa\n', ...
            min(response_avg), max(response_avg));
    else
        response_avg = response_raw;
        fprintf('  警告：无有效内部 DOF，退回单点观测。\n');
    end
else
    response_avg = response_raw;
    fprintf('  警告：无内部节点，退回单点观测。\n');
end

fprintf('正在计算特征频率以便对比…\n');
n_modes = 5;
[eigenvalues, eigenvectors, computed_eigenfreqs] = solveEigenvalueProblem(K, M, air, n_modes);
fprintf('  数值特征频率：');
fprintf('%.1f Hz ', computed_eigenfreqs(1:min(5, length(computed_eigenfreqs))));
fprintf('\n');

fprintf('  各特征频率附近频响幅值：\n');
for i = 1:min(5, length(computed_eigenfreqs))
    [~, closest_idx] = min(abs(frequencies - computed_eigenfreqs(i)));
    response_at_eigen = response_avg(closest_idx);
    fprintf('    %.1f Hz：响应 = %.2e Pa（网格点 %.1f Hz）\n', ...
        computed_eigenfreqs(i), response_at_eigen, frequencies(closest_idx));
end
fprintf('\n');

response_avg_norm = response_avg ./ omega_vec;

peak_locs = [];
resonance_freqs = [];

try
    mean_response = mean(response_avg);
    median_response = median(response_avg);
    threshold = median_response * 1.1;
    [peaks, peak_locs] = findpeaks(response_avg, 'MinPeakHeight', threshold, ...
        'MinPeakProminence', mean_response * 0.05);
    resonance_freqs = frequencies(peak_locs);
catch
    peak_locs = [];
    for i = 2:length(response_avg)-1
        if response_avg(i) > response_avg(i-1) && ...
           response_avg(i) > response_avg(i+1)
            peak_locs = [peak_locs, i];
        end
    end
    resonance_freqs = frequencies(peak_locs);
end

if isempty(peak_locs)
    fprintf('  未找到明显局部峰，在特征频率附近再搜索…\n');
    for i = 1:min(5, length(computed_eigenfreqs))
        if computed_eigenfreqs(i) <= max(frequencies)
            [~, closest_idx] = min(abs(frequencies - computed_eigenfreqs(i)));
            if closest_idx > 1 && closest_idx < length(response_avg)
                if response_avg(closest_idx) > response_avg(closest_idx-1) && ...
                   response_avg(closest_idx) > response_avg(closest_idx+1)
                    peak_locs = [peak_locs, closest_idx];
                    resonance_freqs = [resonance_freqs, frequencies(closest_idx)];
                else
                    mean_nearby = mean([response_avg(max(1, closest_idx-5):min(length(response_avg), closest_idx+5))]);
                    if response_avg(closest_idx) > mean_nearby * 1.2
                        peak_locs = [peak_locs, closest_idx];
                        resonance_freqs = [resonance_freqs, frequencies(closest_idx)];
                    end
                end
            end
        end
    end
end

response = response_avg_norm;
response_raw = response_avg;

fprintf('=== 频响结果摘要 ===\n');
fprintf('共振峰（检出）：\n');
for i = 1:min(5, length(resonance_freqs))
    fprintf('  峰 %d：%.2f Hz\n', i, resonance_freqs(i));
end
fprintf('\n');

expected_freqs = [50, 100, 150, 200, 250];
fprintf('参考纵模频率（3.4 m 硬壁管）：');
fprintf('%.0f Hz ', expected_freqs(1:min(5, length(expected_freqs))));
fprintf('\n\n');

%% 频响图
fig_freq = figure;
set(fig_freq, 'Position', [100, 100, 1200, 600]);

subplot(1, 2, 1);
semilogy(frequencies, response_raw, 'b-', 'LineWidth', 2);
hold on;
for i = 1:min(5, length(expected_freqs))
    xline(expected_freqs(i), 'r--', 'LineWidth', 1.5);
end
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('声压幅值 (Pa)', 'FontSize', 12);
title('频响（原始幅值）', 'FontSize', 14);
legend('响应', '参考纵模', 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

subplot(1, 2, 2);
plot(frequencies, response, 'b-', 'LineWidth', 2);
hold on;
for i = 1:min(5, length(expected_freqs))
    xline(expected_freqs(i), 'r--', 'LineWidth', 1.5);
end
if ~isempty(resonance_freqs) && exist('peak_locs', 'var') && ~isempty(peak_locs)
    plot(resonance_freqs, response(peak_locs), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r');
    legend_str = {'归一化响应', '参考纵模', '检出峰'};
    fprintf('检出 %d 个共振峰：', length(resonance_freqs));
    fprintf('%.1f Hz ', resonance_freqs);
    fprintf('\n');
else
    legend_str = {'归一化响应', '参考纵模'};
    fprintf('未检出明显共振峰（曲线过平或阈值过高）。\n');
end
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('归一化响应 (Pa·s/rad)', 'FontSize', 12);
title('归一化频响（峰标记）', 'FontSize', 14);
legend(legend_str, 'Location', 'best');
grid on;
set(gca, 'FontSize', 11);

%% 谐波声场（选定频率）
fprintf('========================================\n');
fprintf('步骤 3：谐波声场求解与可视化\n');
fprintf('========================================\n');

test_frequencies = [50, 100, 150];
fprintf('在以下频率求解声压场：');
fprintf('%.0f Hz ', test_frequencies);
fprintf('\n');

p_fields = zeros(N.n, length(test_frequencies));

for i = 1:length(test_frequencies)
    freq = test_frequencies(i);
    omega = 2 * pi * freq;
    k = omega / air.c;
    
    f = addAcousticSource(nodes, source_position, source_strength, omega, air);
    
    A = K - (k^2) * M;
    p_fields(:, i) = solveLinearSystem(A, f, 'direct');
    
    fprintf('  第 %d 个频率：%.0f Hz — 已求解\n', i, freq);
end
fprintf('完成。\n\n');

fprintf('正在计算声压级等声学量…\n');
SPL_fields = zeros(N.n, length(test_frequencies));
for i = 1:length(test_frequencies)
    omega = 2 * pi * test_frequencies(i);
    [SPL_fields(:, i), ~, ~] = computeAcousticQuantities(...
        p_fields(:, i), nodes, elements, B, omega, air);
end
fprintf('完成。\n\n');

fprintf('正在绘制声场…\n');
for i = 1:length(test_frequencies)
    fprintf('  频率 %d：%.0f Hz\n', i, test_frequencies(i));
    visualizeAcousticField(nodes, elements, p_fields(:, i), ...
        test_frequencies(i), conf, 1);
    pause(0.5);
end
fprintf('完成。\n\n');

%% 结束
fprintf('========================================\n');
fprintf('频响与声场示例运行结束。\n');
fprintf('========================================\n');
fprintf('输出：频响曲线、各频率声场与 SPL 分布。\n');
fprintf('\n');

function node_idx = findClosestNode(nodes, position)
    node_coords = nodes.all(:, 1:2);
    distances = sqrt(sum((node_coords - repmat(position(1:2), size(node_coords, 1), 1)).^2, 2));
    [~, node_idx] = min(distances);
end
