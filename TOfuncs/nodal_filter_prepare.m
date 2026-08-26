function [node_weight_sum, node_weight_kernel, weight_elem2node] = nodal_filter_prepare(nelx, nely, rmin, passive)
% 生成基于结点的过滤器矩阵及其行和
% 输入参数：
%   nelx : 设计域在x方向的单元数量
%   nely : 设计域在y方向的单元数量
%   rmin : 过滤器半径
%   passive: 设计域标记 (可选，等于1为非设计域，默认值全0)
% 输出参数：
%   node_weight_sum    : 结点的权重和矩阵，大小为(nely x nelx)
%   node_weight_kernel : 结点的权重核矩阵，大小为(ceil(rmin - 0.5) * 2)^2
%   weight_elem2node   : 单元值到节点值的转换矩阵，大小为(num_nodes x num_elems)
% 用法：
% node_alpha_i = sen_weight_elem2node * alpha_i(:); % 灵敏度：单元 -> 结点
% node_alpha_i = conv2(reshape(node_alpha_i, nely+1, nelx+1), sen_node_weight_kernel, "same"); % 算分子：基于结点进行求和操作
% alpha_i = node_alpha_i(1:end-1, 1:end-1) ./ sen_node_weight_sum; % 除以分母：过滤操作完成，因为是偶数核，核的中心认为是中间偏左上，所以最下一行和最后一列是多余的

% 最后修改: 
% H Yingvar 2026-01-07

%% 初始化和参数处理
if nargin < 4
    passive = false(nely, nelx);
else
    passive = reshape(passive, nely, nelx);
    passive = logical(passive);
end

num_elems = nelx * nely;
num_nodes = (nelx + 1) * (nely + 1);
node_mat = reshape(1:num_nodes, nely + 1, nelx + 1); % 节点编号矩阵
% 左下节点
n1 = node_mat(2:nely+1, 1:nelx);
% 右下节点
n2 = node_mat(2:nely+1, 2:nelx+1);
% 右上节点
n3 = node_mat(1:nely, 2:nelx+1);
% 左上节点
n4 = node_mat(1:nely, 1:nelx);
elem_mat = [n1(:) n2(:) n3(:) n4(:)]; % 第i行表示单元i对应的结点编号，左下开始逆时针转
% in_design_elems_index = find(passive == 0); % 设计域内单元的索引
% in_design_elem_mat = elem_mat(~passive,:);
% in_design_nodes = unique(in_design_elem_mat);
% non_design_nodes = setdiff((1:num_nodes)', in_design_elem_mat); % 非设计域
node_passive = conv2(~passive, [1 1;1 1]) == 0; % 标记结点是否在设计域内的矩阵

%% 1. 建立单元到节点的平均矩阵
weight_elem2node = sparse(elem_mat(:), repmat((1:num_elems)', 4, 1), 1, num_nodes, num_elems);
weight_elem2node(node_passive, :) = 0;
weight_elem2node(:, passive) = 0;
sum_divide_by = sum(weight_elem2node, 2);
sum_divide_by(sum_divide_by == 0) = 1;
weight_elem2node = weight_elem2node ./ sum_divide_by;

%% 2. 建立节点的卷积核，kernel_r的含义是（1/2 * 矩阵尺寸）
kernel_r = ceil(rmin-0.5) + 1/2 * (-1 + sign(rmin - hypot(ceil(rmin-0.5)-0.5, 0.5))); % 先估计一个保守值，然后判断离半径最近的结点是否在过滤域内，若在则取该值；若不在则减1
[dx, dy] = meshgrid(.5 - kernel_r : kernel_r-0.5);
node_weight_kernel = max(0, rmin - hypot(dx, dy));

%% 3. 创建最终的权重和矩阵
node_weight_sum = conv2(~node_passive, node_weight_kernel, "same");
node_weight_sum = node_weight_sum(1:end-1, 1:end-1);
node_weight_sum(passive) = 1e9; % 1. 防止除以0时异常，2. 使非设计域的灵敏度趋近于0

end
