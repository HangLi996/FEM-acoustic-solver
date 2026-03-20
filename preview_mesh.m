% preview_mesh.m — 快速预览网格与边界
%
% 从 Gmsh .msh 读入网格，识别边界节点，并调用 plotMesh 绘制。
% 适合在运行完整分析前检查几何与剖分是否合理。
%
% 用法：在仓库根目录执行 preview_mesh

clear all; close all; clc;
addpath('functions');
addpath('mesh');
config;

%% 用户可改：网格路径
file.mesh = 'mesh/duct/duct4l.msh';

%% 读入网格与边界
[nodes, elements, dim, N] = msh2mat(file);
[B, N] = getBoundary(file, conf, N, nodes, elements);

%% 绘图
fig_mesh = figure;
movegui(fig_mesh, 'northwest');
hold all;
plotMesh(fig_mesh, nodes, elements, conf, B, [1 1 1], 'k', [0 90], 'northwest');
title('网格与边界预览');
