%% testDK.m — 排序函数（带宽/填充）对比测试
% 对同一个缩减刚度矩阵 Kr = K(free,free)，分别用不同排序重排后再做稀疏
% Cholesky，对比：排序耗时、分解耗时、分解因子非零数 nnz(L)、填充比、带宽。
% 并校验各排序下求得的解一致（验证置换回代正确）。
%
% 网格规模说明：3D 列主序编号下 z 方向带宽巨大，natural（不重排）chol 的填充
% 约为 N * 3*(nelx+1)*(nely+1)，随网格急剧膨胀：
%   60x60x20（23 万自由度）→ 填充 ~37.3GB → 本机 16GB 内存直接报错
%   100x100x20（64 万自由度）→ 填充估算 ~157GB，且组集就要 ~4.6GB，无法运行
% 为能安全测出各排序的差异（含会爆的 natural），本脚本默认用 30x30x20
% （natural 填充 ~2.8GB，仍可跑完）。想测更大规模可上调 nelx/nely，
% 但 natural 会先内存超限，此时需注释掉 natural 那一项。

clear; clc; close all;

%% 1. 网格与材料参数（与 pnormBesoSubD_3D / test.m 保持一致）
nelx = 30;  % x 方向单元数（30x30x20 才安全；100 会内存超限，见头部说明）
nely = 30;  % y 方向单元数
nelz = 10;  % z 方向单元数

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

%% 4. 排序函数对比：natural / symrcm / symamd / amd / dissect
pdof = BC(:, 1);
free = setdiff((1:ndof)', pdof);
Kr   = K(free, free);
n    = numel(free);
b0   = F(free);                       % 齐次 BC 的右端项
nnzK = nnz(Kr);

orderings = {'natural', 'symrcm', 'symamd', 'amd', 'dissect'};
fprintf('\n==== 排序函数对比（自由度=%d, Kr nnz=%d）====\n', n, nnzK);
fprintf('%-9s %10s %10s %13s %9s %10s %12s\n', ...
        '排序', '排序耗时s', '分解耗时s', 'L nnz', '填充比', '带宽', '解误差');
dref = [];
for i = 1:numel(orderings)
    name = orderings{i};
    switch name
        case 'natural', t0 = tic; p = (1:n)';          t_order = toc(t0); % 不重排（基线）
        case 'symrcm',  t0 = tic; p = symrcm(Kr);      t_order = toc(t0); % 带宽最小化
        case 'symamd',  t0 = tic; p = symamd(Kr);      t_order = toc(t0); % 填充最小化
        case 'amd',     t0 = tic; p = amd(Kr);         t_order = toc(t0); % 填充最小化
        case 'dissect', t0 = tic; p = dissect(Kr);     t_order = toc(t0); % 嵌套剖分填充最小化
    end
    t0 = tic;
    L  = chol(Kr(p, p), 'lower');
    t_chol = toc(t0);
    d = zeros(n, 1);
    d(p) = L' \ (L \ b0(p));         % 置换空间回代
    if isempty(dref), dref = d; err = 0; else, err = norm(d - dref, inf); end
    bw = bandwidth(Kr(p, p), 'lower');
    fprintf('%-9s %10.3f %10.3f %13d %9.1f %10d %12.2e\n', ...
            name, t_order, t_chol, nnz(L), nnz(L)/nnzK, bw, err);
end
fprintf('解误差：各排序相对 natural 解的无穷范数差异，理论应全为 0\n');
fprintf('带宽：重排后矩阵 Kr(p,p) 的下带宽（非零元离主对角线的最远距离）\n');
fprintf('填充比：nnz(L)/nnz(Kr)，越小代表分解产生的新非零元越少、越省内存\n');