function plot_metrics(metric_names, metric_data, save_path, title_str)
% plot_metrics 绘制多个 metric 随迭代次数变化的曲线并保存数据与图片
% 可选参数:
%   save_path: 保存路径，默认 `pwd`
%   title_str: 绘图标题，默认使用保存目录的上级目录名
%
% 用法:
%   plot_metrics(metric_names, metric_data, save_path)
% 参数:
%   metric_names: cell array of strings, 每一列对应一个 metric 名称
%   metric_data: 数值矩阵, 行表示迭代次数, 列对应 metric
%   save_path: 保存路径，可以是目录或文件路径（函数会在目录下保存文件）
%
% 输出:
%   在 save_path 指定的目录下保存:
%     - metrics_plot.png: 绘图图片
%     - metrics_plot.fig: MATLAB 图形文件
%     - metrics_data.csv: 包含表头（metric_names）和数值数据的 CSV

if nargin < 2
    error('需要至少 2 个输入参数：metric_names, metric_data');
end

% 处理可选参数
if nargin < 3 || isempty(save_path)
    save_path = pwd;
end


% 检查 metric_names
if ischar(metric_names) || isstring(metric_names)
    % 支持 char 向量/矩阵 与 string 标量/数组
    metric_names = cellstr(metric_names);
end
if ~iscell(metric_names)
    error('metric_names 应为 cell array of strings');
end

% 检查 metric_data
if ~isnumeric(metric_data) || ~ismatrix(metric_data)
    error('metric_data 应为二维数值矩阵');
end

[nIter, nMetric] = size(metric_data);
if numel(metric_names) ~= nMetric
    error('metric_names 的数量必须等于 metric_data 的列数');
end

% 处理保存目录
if exist(save_path, 'dir') == 7
    outdir = save_path;
else
    [p, ~, ~] = fileparts(save_path);
    if isempty(p)
        outdir = pwd;
    else
        outdir = p;
    end
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end


% 准备颜色和标记
colors = lines(nMetric);
markers = {'o','s','d','^','v','>','<','p','h','+','x','*','.'};

% 创建图形并为每个 metric 使用子图
fig = figure('Color','w');

% 子图布局：近似正方形
nrows = ceil(sqrt(nMetric));
ncols = ceil(nMetric / nrows);
% 自适应标记策略：每条曲线上最多显示的标记数量（可根据需要调整）
maxMarkers = 30;
for i = 1:nMetric
    ax = subplot(nrows, ncols, i);
    mk = markers{mod(i-1, numel(markers)) + 1};
    % 先画连线
    plot(1:nIter, metric_data(:,i), '-','Color', colors(i,:), 'LineWidth', 1.5);
    hold(ax, 'on');
    % 根据数据点数量计算标记间隔，保证每条曲线的标记数不超过 maxMarkers
    markerStep = max(1, ceil(nIter / maxMarkers));
    markerIdx = 1:markerStep:nIter;
    % 单独绘制标记（只在稀疏索引处绘制），避免点过密影响观感
    plot(markerIdx, metric_data(markerIdx, i), 'LineStyle', 'none', ...
        'Marker', mk, 'Color', colors(i,:), 'MarkerSize', 6, 'MarkerFaceColor', colors(i,:));
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'Iteration', 'Interpreter', 'latex');
    ylabel(ax, metric_names{i}, 'Interpreter', 'latex');
    set(ax, 'Box', 'on');
end

if nargin > 3
   sgtitle(title_str, 'FontSize', 16, 'FontName', 'Times New Roman', 'Fontweight', 'Bold');
end
% 保存图片（不保存 .fig）
pngfile = fullfile(outdir, 'metrics_plot.png');
try
    saveas(fig, pngfile);
catch
    warning('保存图片失败。');
end

% 保存 CSV: 先写表头，再写数值
csvfile = fullfile(outdir, 'metrics_data.csv');
fid = fopen(csvfile, 'w');
if fid == -1
    warning('无法创建 CSV 文件：%s', csvfile);
    return;
end
% 写表头，必要时对包含逗号或双引号的 header 添加双引号
for i = 1:nMetric
    header = metric_names{i};
    header = strrep(header, '"', '""');
    if contains(header, ',') || contains(header, '"')
        fprintf(fid, '"%s"', header);
    else
        fprintf(fid, '%s', header);
    end
    if i < nMetric
        fprintf(fid, ',');
    else
        fprintf(fid, '\n');
    end
end

% 写数据行
fmt = repmat('%.15g,', 1, nMetric);
fmt(end) = newline; % 最后一个逗号替换为换行 (使用单字符换行符)
for r = 1:nIter
    line = metric_data(r, :);
    fprintf(fid, fmt, line);
end
fclose(fid);

% 输出保存信息
fprintf('Saved plot: %s\n', pngfile);
fprintf('Saved data: %s\n', csvfile);
end
