%% benchmark_solveq_vs_solveq2.m
% 对比 solveq 与 solveq2 求解「带 BC 线性系统」的速度
%   连续两次求解 K*d = f（相同的 K、F、BC），分别计时：
%     - solveq :  第1次完整分解（Cholesky+dissect 排序）并缓存；第2次复用缓存，仅三角回代（O(1)）
%     - solveq2:  每次都完整分解（\），无缓存
% 两组求解使用完全相同的 K、F、BC，并在最后校验结果一致性。
%
% 网格规模说明：3D 列主序编号下，直接 chol（solveq 旧版）与 mldivide（solveq2）
% 在 z 方向填充灾难，60x60x20（约 23 万自由度）在本机（16GB RAM）会内存超限
% （报错 232848x232848 / 37.3GB）。solveq 已改用 dissect 嵌套剖分排序解决；
% 但 solveq2 的 "\" 仍无法处理 60x60x20，故本脚本默认用 40x40x20（约 10.6 万
% 自由度），两者都能在数分钟内跑完。若只测 solveq，可把 nelx/nely 调到 60。

clear; clc; close all;

%% 1. 网格与材料参数（与 pnormBesoSubD_3D / test.m 保持一致）
nelx = 10;  % x 方向单元数
nely = 10;  % y 方向单元数
nelz = 20;  % z 方向单元数
subD_size = 0.5;

E0  = 2.1e5; % 杨氏模量
nu  = 0.3;   % 泊松比
ex  = [0, 1, 1, 0, 0, 1, 1, 0]; % 单元节点 x 坐标
ey  = [0, 0, 1, 1, 0, 0, 1, 1]; % 单元节点 y 坐标
ez  = [0, 0, 0, 0, 1, 1, 1, 1]; % 单元节点 z 坐标
ep  = 2;                        % 2x2x2 高斯积分点
D0  = hooke(4, E0, nu);         % 3D 材料本构矩阵
k0  = soli8e(ex, ey, ez, ep, D0); % 单元刚度矩阵

%% 2. 设计域（长方体）、载荷、边界条件
% 全实心长方体：左端面（x=1）全部自由度固定，右端面（x=nelx）施加 x 方向均匀拉力
x = ones(nely, nelx, nelz) > 0;   % 设计变量：全 1（实心）

nodeMat = reshape(1:(1+nelx)*(1+nely)*(1+nelz), nely+1, nelx+1, nelz+1);

% 左端面固定：nodeMat(:,1,:) 为 x=1 面上全部节点，固定其 3 个平动自由度
fixed_nodes = unique(nodeMat(:, 1, :));
fixed_dofs  = kron(3*fixed_nodes, ones(3,1)) + repmat([-2 -1 0]', numel(fixed_nodes), 1);
BC = [fixed_dofs, zeros(numel(fixed_dofs), 1)];

% 右端面施加 x 方向均匀拉力：总力 1e3，平均分到 x=nelx 面各节点
load_nodes = unique(nodeMat(:, end, :));
load_dofs  = 3*load_nodes - 2;    % 各节点的 x 平动自由度
F = sparse(load_dofs, 1, 1e3/numel(load_dofs), 3*(1+nelx)*(1+nely)*(1+nelz), 1);

%% 3. 组集全局刚度矩阵 K（流程与 structure_FEA_subD_3D 相同）
[nely, nelx, nelz] = size(x);
num_elems = nelx * nely * nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);

nodeMat  = reshape(1:(1+nelx)*(1+nely)*(nelz+1), nely+1, nelx+1, nelz+1);
eNodeVec = reshape(nodeMat(2:end, 1:end-1, 2:end), num_elems, 1);
eNodeMat = repmat(eNodeVec,1,8) + repmat([0 nely+[1 0] -1 -1*(nely+1)*(nelx+1)+[0 nely+[1 0] -1]], num_elems, 1);
edofMat  = kron(eNodeMat, ones(1,3)) * 3 + repmat([-2 -1 0], num_elems, 8);

iK = reshape(kron(edofMat, ones(24,1))', 24*24*num_elems, 1);
jK = reshape(kron(edofMat, ones(1,24))', 24*24*num_elems, 1);
sK = reshape(k0(:)*max(1e-9, x(:))', 24*24*num_elems, 1);
K  = sparse(iK, jK, sK);
clear nodeMat eNodeVec eNodeMat edofMat iK jK sK; % 释放组集临时矩阵

fprintf('问题规模：%d 个单元，%d 个自由度，%d 条 BC\n', num_elems, ndof, size(BC,1));

%% 4. solveq：连续两次求解（第1次分解+缓存，第2次复用缓存）
clear solveq; % 清空可能残留的持久缓存，保证公平
t0 = tic;
[d1a, Q1a] = solveq(K, F, BC);        % 第 1 次：分解 + 三角回代，并缓存
t_solveq_1 = toc(t0);

t0 = tic;
[d1b, Q1b] = solveq(K, F, BC, true);  % 第 2 次：复用缓存，仅三角回代
t_solveq_2 = toc(t0);

% %% 5. solveq2：连续两次求解（每次都完整分解，无缓存）
% t0 = tic;
% [d2a, Q2a] = solveq2(K, F, BC);       % 第 1 次
% t_solveq2_1 = toc(t0);

% t0 = tic;
% [d2b, Q2b] = solveq2(K, F, BC);       % 第 2 次
% t_solveq2_2 = toc(t0);

% %% 6. 结果一致性校验（防止对比失真）
% tol = 1e-8;
% err_d = max(norm(d1a - d2a, inf), norm(d1b - d2b, inf));
% err_Q = max(norm(Q1a - Q2a, inf), norm(Q1b - Q2b, inf));
% if err_d < tol && err_Q < tol
%     fprintf('✔ 结果一致：d 最大误差 %.2e，Q 最大误差 %.2e\n', err_d, err_Q);
% else
%     warning('结果不一致！d 最大误差 %.2e，Q 最大误差 %.2e', err_d, err_Q);
% end

%% 7. 汇总输出
fprintf('\n========== 单次调用耗时（秒） ==========\n');
fprintf('  solveq  第1次（分解+回代+缓存） : %8.4f s\n', t_solveq_1);
fprintf('  solveq  第2次（复用缓存回代）   : %8.4f s\n', t_solveq_2);
% fprintf('  solveq2 第1次（完整分解）       : %8.4f s\n', t_solveq2_1);
% fprintf('  solveq2 第2次（完整分解）       : %8.4f s\n', t_solveq2_2);

fprintf('\n========== 连续两次总耗时（秒） ==========\n');
fprintf('  solveq  : %8.4f s\n', t_solveq_1 + t_solveq_2);
% fprintf('  solveq2 : %8.4f s\n', t_solveq2_1 + t_solveq2_2);

fprintf('\n========== 加速比 ==========\n');
% fprintf('  第2次调用   solveq2/solveq : %6.2f x\n', t_solveq2_2 / t_solveq_2);
% fprintf('  连续两次     solveq2/solveq : %6.2f x\n', (t_solveq2_1 + t_solveq2_2) / (t_solveq_1 + t_solveq_2));
