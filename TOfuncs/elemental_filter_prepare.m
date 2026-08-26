function [Hs, h] = elemental_filter_prepare(nelx, nely, nelz, rmin, passive)
% [Hs, h] = elemental_filter_prepare(nelx, nely, nelz, rmin, passive)
% 生成基于单元的过滤器矩阵H及其行和Hs（支持二维与三维）
% 输入参数：
%   nelx : 设计域在x方向的单元数量
%   nely : 设计域在y方向的单元数量
%   nelz : 设计域在z方向的单元数量（二维问题时输入1或空[]）
%   rmin : 过滤器半径
%   passive: 设计域标记 (可选，等于1为非设计域，默认值全0；其维度决定二维或三维处理)
% 输出参数：
%   h  : 过滤器卷积核矩阵，二维时大小为(2 * (ceil(rmin)-1) +1)^2，三维时为
%        (2 * (ceil(rmin)-1) +1)^3，对应以某单元为中心rmin范围内的所有单元对应的权重
%   Hs : 过滤器权重和矩阵，大小与设计域一致（二维为 nely x nelx，三维为 nely x nelx x nelz），
%        以单元(ely,elx[,elz])为中心对应的rmin范围内权重总和
% 用法：
% 二维: alpha_i = conv2(alpha_i, h, 'same') ./ Hs;
% 三维: alpha_i = convn(alpha_i, h, 'same') ./ Hs;

% 最后修改:
% H Yingvar 2026-08-12
% change log:
% 1. 支持三维输入（convn 与三维过滤核）
% 2. 参数顺序调整为 (nelx, nely, nelz, rmin, passive)，二维问题时 nelz 输入 1 或空[]

%% 初始化和参数处理
% nelz 为空时按二维处理（rmin 为必选参数，因此 nelz 必定被传入）
if isempty(nelz)
    nelz = 1;
end

if nargin < 5 || isempty(passive)
    % passive未提供：依据nelz决定维度
    if nelz > 1
        passive = false(nely, nelx, nelz);
    else
        passive = false(nely, nelx);
    end
else
    % passive 已按设计域形状（二维或三维）输入，无需 reshape
    passive = logical(passive);
end

is3D = ndims(passive) == 3; % 是否为三维设计域

%% 计算过滤核，kernel_r的含义是（1/2 * 矩阵尺寸 -1）
kernel_r = ceil(rmin);
if is3D
    v = 1-kernel_r : kernel_r-1;
    [dx, dy, dz] = ndgrid(v, v, v);
    h = max(0, rmin - sqrt(dx.^2 + dy.^2 + dz.^2)); % 三维欧氏距离权重
else
    [dx, dy] = meshgrid(1-kernel_r:kernel_r-1);
    h = max(0, rmin - hypot(dx, dy));  % 更高效的距离计算
end

%% 创建最终的权重和矩阵
if is3D
    Hs = convn(~passive, h, 'same'); % 三维卷积
else
    Hs = conv2(~passive, h, "same");
end
Hs(logical(passive)) = 1e9; % 1. 防止除以0时异常，2. 使非设计域的灵敏度趋近于0

end

