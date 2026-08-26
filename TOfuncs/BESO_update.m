function x_new = BESO_update(x_old, alpha_i, vol_new, passive, AR_max, is_soft_kill)
% 功能：
% 按照BESO算法的规则根据灵敏度和目标体积更新设计变量
% 语法：
% x_new = BESO_update(x_old, alpha_i, vol_new, AR_max) 
% 输入：
%   x_old - 当前设计变量 (去除标记为passive的元素后输入）
%   alpha_i - 灵敏度 (去除标记为passive的元素后输入）
%   vol_new - 目标体积
%   passive - 设计域标记 (可选，等于1为非设计域，默认值全0)
%   AR_max - 最大准加比 (可选，默认值0.02)
%   is_soft_kill - 确定是否为soft-kill，默认为false（hard-kill)
% 输出：
%   x_new - 更新后的设计变量

% 最后修改: 
% H Yingvar 2026-02-23

%% 初始化参数
if nargin < 4
    passive = zeros(size(x_old)) > 0; % 默认没有非设计域
    AR_max = 0.2; % 默认最大准加比
    is_soft_kill = false;
elseif nargin < 5
    AR_max = 0.2; % 默认最大准加比
    is_soft_kill = false;
elseif nargin < 6
    is_soft_kill = false; % 默认为false（hard-kill)
end

if is_soft_kill
    x_old = x_old > 0.99;
end

passive = logical(passive);
x_old_process = x_old(~passive);
alpha_process = alpha_i(~passive); % 去除非设计域后的灵敏度
% vol_domain = length(x_old_process); % 设计域总体积
vol_old = nnz(x_old_process); % 当前设计的体积
x_new_process = zeros(size(x_old_process)); % 初始化更新后的设计变量
x_new = x_old; % 初始化更新后的设计变量

%% 仅处理设计域中的单元
[~, alpha_sort_indicas] = sort(alpha_process, 'descend');
x_new_process(alpha_sort_indicas(1:round(vol_new))) = 1; % 高于阈值的部分设为1，满足体积约束
vol_add = nnz(x_new_process - x_old_process == 1);
AR = vol_add / vol_old;
if AR > AR_max
    vol_add = AR_max * vol_old;
    alpha_sort_indicas_void = alpha_sort_indicas(x_old_process(alpha_sort_indicas) == 0); % 对空洞单元的灵敏度排序索引
    x_new_process = x_old_process; % 重置更新后的设计变量
    x_new_process(alpha_sort_indicas_void(1:round(vol_add))) = 1; % 调整添加体积以满足最大准加比
    elem_remove_count = round(vol_old - vol_new + vol_add);
    alpha_sort_indicas_solid = alpha_sort_indicas(x_old_process(alpha_sort_indicas) == 1); % 对实心单元的灵敏度排序索引
    x_new_process(alpha_sort_indicas_solid(end - elem_remove_count + 1:end)) = 0; % 调整删除体积以满足体积约束
end
x_new(~passive) = x_new_process; % 将更新后的设计变量放回原始设计变量矩阵中
if nnz(x_new_process) - round(vol_new) > 1
    disp("体积超了，检查BESO_update")
end
if is_soft_kill
    x_new = double(x_new);
    x_new(x_new == 0 & ~passive) = 1e-9;
end

end
