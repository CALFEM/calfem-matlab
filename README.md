# TOptLab · 拓扑优化实验室

> **命名说明**：`TOptLab` = **T**opology **Opt**imization **Lab**oratory（拓扑优化实验室），直指核心用途——算法快速验证；其中 `TOpt` 是拓扑优化（Topology Optimization）领域的通用缩写，易检索、好传播。

## 项目介绍

本项目旨在构建一个**适用于拓扑优化算法快速验证的功能函数库**，仅使用纯 MATLAB 语法实现，不依赖任何第三方工具箱。当前规划包含两个核心模块：

- **有限元运算函数库** —— 见文件夹 `fem`
- **拓扑优化函数库** —— 见文件夹 `TOfuncs`

## 目录结构

| 目录       | 说明               |
| ---------- | ------------------ |
| `fem`      | 有限元运算函数库   |
| `TOfuncs`  | 拓扑优化函数库     |
| `test`     | 测试脚本           |

## 有限元模块发展情况

基础文件拷贝自 [CALFEM](https://calfem.com)，一个面向 MATLAB 的有限元工具箱。

> A finite element toolbox for MATLAB. © Division of Structural Mechanics and Division of Solid Mechanics, Lund University

### 更新日志

- **`assem()`**：提高了大型刚度矩阵的组装速度，输入参数有修改。
- **`extract_ed()`**：提高了从全局位移向量中提取单元位移的速度，输入参数有修改。
- **`solveq()`**：优化了超大型平衡方程求解的内存使用与运算速度。
- **`soli8s()`**：新增 `ir=1`（1×1×1 中心点）积分规则；修正切应力分量顺序为 `[sigxy sigyz sigxz]`；移除积分点循环，改为一次性向量化计算所有积分点；当所有单元共享同一几何时完全无循环，显著提升速度。

## 拓扑优化模块发展情况
目前所有函数均由作者本人编写，模块仍在持续建设中……

### 更新日志

- **`BESO_update()`**：按照 BESO 算法规则，基于灵敏度与目标体积更新设计变量；支持 hard-kill / soft-kill 两种模式，并可通过 `passive` 标记非设计域。
- **`elemental_filter_prepare()`**：生成基于单元的过滤器卷积核 `h` 及其权重和矩阵 `Hs`，支持二维与三维设计域（分别使用 `conv2` / `convn`）。
- **`nodal_filter_prepare()`**：生成基于节点的过滤器矩阵，包含节点权重和、权重核，以及单元到节点的值转换矩阵。
- **`plot_metrics()`**：绘制多个指标随迭代次数变化的曲线，并自动保存图片与数据文件（PNG、FIG、CSV）。
- **`display_2D()`**：显示二维拓扑标量场（第一个参数为密度场 `rho`），名称-值对接口与 `display_3D` 一致；`rho <= 0.5` 的孔洞单元自动隐藏，支持指定着色场、掩膜、颜色范围、色图与 colorbar，并可选仿 FPTO2D 的光滑边界显示。
- **`display_3D()`**：显示三维拓扑，由 `Smooth` 开关切换两种模式——体素方块按标量场逐单元着色，或 `isosurface` 平滑显示（top3D125 蓝色主体配色）；支持指定着色场、掩膜、颜色范围、色图与 colorbar。
- **`model2STL()`**：由模型生成 STL 文件（立方体 / 平滑 isosurface 双模式），参考 Top3D 扩展 `Top3dSTL_v3.m`，用于三维打印或后处理；支持二进制/ASCII 格式、密度阈值、面片颜色与透明度、单元边长等选项，并可选择是否同时可视化。
