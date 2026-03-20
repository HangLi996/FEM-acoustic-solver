% run_eigen_modes.m — 声学有限元：模态分析示例
%
% 本脚本演示 FEM-acoustic-solver 的典型工作流：
%   1. 读取 Gmsh 网格，计算单元雅可比与行列式（单元质量检查）
%   2. 组装全局刚度矩阵 K 与质量矩阵 M
%   3. 求解广义特征值问题并可视化模态
%
% 用法：在仓库根目录执行 run_eigen_modes
%
% 依赖：functions/ 下工具函数；默认网格见下方 file.mesh

clear all; close all; clc;
addpath('functions');
addpath('mesh');
config;

%% 用户可改：网格文件（相对仓库根目录）
file.mesh = 'mesh/duct/duct4l.msh';

fprintf('========================================\n');
fprintf('  FEM 声学求解器 — 模态分析示例\n');
fprintf('========================================\n');

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

%% 步骤一：雅可比矩阵与行列式（单元映射健康检查）
fprintf('========================================\n');
fprintf('步骤 1：各单元雅可比与行列式\n');
fprintf('========================================\n');

Jacobian = cell(N.el, 1);
detJacobian = zeros(N.el, 1);

fprintf('正在计算各单元雅可比…\n');
for el = 1:N.el
    el_nodes = elements.all(el, :);
    el_nodes = el_nodes(el_nodes ~= 0);
    node_coords = nodes.all(nodes.nind(el_nodes, 2), 1:2);
    
    element_type = elements.type(el);
    
    % 在单元中心处计算雅可比
    if element_type == 1 || element_type == 8
        [J, detJ] = computeJacobian(0, [], node_coords, element_type);
    else
        [J, detJ] = computeJacobian(0, 0, node_coords, element_type);
    end
    
    Jacobian{el} = J;
    detJacobian(el) = detJ;
    
    if detJ < 0
        warning('单元 %d 雅可比行列式为负：%.6f', el, detJ);
    end
end

fprintf('完成。\n');
fprintf('  行列式范围：[%.6e, %.6e]\n', ...
    min(detJacobian), max(detJacobian));
fprintf('  负行列式单元数：%d\n', sum(detJacobian < 0));
fprintf('  近似零行列式单元数：%d\n', sum(abs(detJacobian) < 1e-10));
fprintf('\n');

%% 步骤二：组装全局 K、M
fprintf('======================================\n');
fprintf('步骤 2：组装全局刚度矩阵 K 与质量矩阵 M\n');
fprintf('========================================\n');

fprintf('正在计算单元矩阵并组装…\n');
tic;
[K, M] = assembleGlobalMatrices(nodes, elements, N);
time_assembly = toc;

fprintf('完成（耗时 %.4f 秒）\n', time_assembly);
fprintf('  全局刚度矩阵 K：[%d × %d]\n', size(K, 1), size(K, 2));
fprintf('  全局质量矩阵 M：[%d × %d]\n', size(M, 1), size(M, 2));
fprintf('  K 非零元：%d（占满阵 %.2f%%）\n', ...
    nnz(K), 100*nnz(K)/numel(K));
fprintf('  M 非零元：%d（占满阵 %.2f%%）\n', ...
    nnz(M), 100*nnz(M)/numel(M));
fprintf('\n');

% 稀疏结构可视化
fig_sparse = figure;
set(fig_sparse, 'Position', [100, 100, 1200, 500]);
subplot(1, 2, 1);
spy(K);
title('刚度矩阵 K 的稀疏结构');
xlabel('列索引');
ylabel('行索引');

subplot(1, 2, 2);
spy(M);
title('质量矩阵 M 的稀疏结构');
xlabel('列索引');
ylabel('行索引');

%% 步骤三：广义特征值问题 K·p = k²·M·p
fprintf('========================================\n');
fprintf('步骤 3：求解广义特征值问题\n');
fprintf('======================================\n');

fprintf('正在求解 K·p = k²·M·p …\n');
n_modes = 10;
[eigenvalues, eigenvectors, eigenfrequencies] = ...
    solveEigenvalueProblem(K, M, air, n_modes);

fprintf('完成。\n');
fprintf('\n');
fprintf('=== 特征值结果 ===\n');
fprintf('%-8s  %-18s  %-15s  %-15s\n', '模态', '特征值(k²)', '波数 k', '频率(Hz)');
fprintf('%s\n', repmat('-', 1, 62));

for i = 1:min(10, length(eigenfrequencies))
    k = sqrt(eigenvalues(i));
    fprintf('%-8d  %-18.6e  %-15.6f  %-15.2f\n', ...
        i, eigenvalues(i), k, eigenfrequencies(i));
end
fprintf('\n');

% 与一维硬壁管道纵模解析频率对照（长度约 3.4 m 时）
fprintf('=== 频率核对（参考：3.4 m 硬壁管道纵模）===\n');
fprintf('理论纵模约：50 Hz, 100 Hz, 150 Hz, …\n');
fprintf('本算例前 5 个数值频率：');
fprintf('%.2f Hz ', eigenfrequencies(1:min(5, length(eigenfrequencies))));
fprintf('\n\n');

%% 模态振型
fprintf('正在绘制模态振型…\n');
mode_indices = 1:min(6, length(eigenfrequencies));
visualizeModeShapes(nodes, elements, eigenvectors, mode_indices, ...
    eigenfrequencies, conf);
fprintf('完成。\n\n');

%% 频谱柱状图
fig_spectrum = figure;
set(fig_spectrum, 'Position', [200, 200, 800, 600]);
bar(eigenfrequencies(1:min(10, length(eigenfrequencies))));
xlabel('模态序号', 'FontSize', 12);
ylabel('特征频率 (Hz)', 'FontSize', 12);
title('特征频率谱', 'FontSize', 14);
grid on;
set(gca, 'FontSize', 11);

%% 结束
fprintf('========================================\n');
fprintf('模态分析示例运行结束。\n');
fprintf('========================================\n');
fprintf('输出：稀疏结构图、模态振型图、特征频率谱。\n');
fprintf('\n');
