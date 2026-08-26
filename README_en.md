# TOptLab · Topology Optimization Laboratory

> **About the name**: `TOptLab` stands for **T**opology **Opt**imization **Lab**oratory, pointing directly to its core purpose — the rapid validation of topology optimization algorithms. `TOpt` is a widely used abbreviation in the field of topology optimization, making the name easy to search for and easy to share.

## Project Overview

This project aims to build a **function library for the rapid validation of topology optimization algorithms**, implemented using pure MATLAB syntax with no dependency on third-party toolboxes. The current plan consists of two core modules:

- **Finite element library** — see folder `fem`
- **Topology optimization library** — see folder `TOfuncs`

## Directory Structure

| Directory  | Description                        |
| ---------- | ---------------------------------- |
| `fem`      | Finite element library             |
| `TOfuncs`  | Topology optimization library      |
| `test`     | Test scripts                       |

## Finite Element Module

The base files are taken from [CALFEM](https://calfem.com), a finite element toolbox for MATLAB.

> A finite element toolbox for MATLAB. © Division of Structural Mechanics and Division of Solid Mechanics, Lund University

### Changelog

- **`assem()`**: Improved the speed of assembling large stiffness matrices; function signature modified.
- **`extract_ed()`**: Improved the speed of extracting element displacements from the global displacement vector; function signature modified.
- **`solveq()`**: Optimized memory usage and computational speed for solving very large equilibrium systems.

## Topology Optimization Module

All functions are currently written by the author, and the module is still under active development.

### Changelog

- **`BESO_update()`**: Updates design variables according to the BESO algorithm rules based on sensitivities and a target volume; supports both hard-kill and soft-kill schemes, and allows marking non-design domains via `passive`.
- **`elemental_filter_prepare()`**: Builds the element-based filter convolution kernel `h` and its weight-sum matrix `Hs`, supporting 2D and 3D design domains (via `conv2` / `convn` respectively).
- **`nodal_filter_prepare()`**: Builds a node-based filter, including the node weight sums, the weight kernel, and an element-to-node value transformation matrix.
- **`plot_metrics()`**: Plots multiple metrics against iteration number and automatically saves the figures and data files (PNG, FIG, CSV).
