  function [d,Q]=solveq(K,f,bc,reuse)
% a=solveq(K,f)
% [a,Q]=solveq(K,f,bc)
%-------------------------------------------------------------
% PURPOSE
%  Solve static FE-equations considering boundary conditions.
%
% INPUT: K : global stiffness matrix, dim(K)= nd x nd
%        f : global load vector, dim(f)= nd x 1
%        bc : boundary condition matrix
%            dim(bc)= nbc x 2, nbc : number of b.c.'s
%        reuse : (optional) 若为 true：复用缓存中的分解（O(1)，仅三角回代）；
%            调用方必须保证 K 与 bc 与最近一次 solveq 调用完全一致（详见下方说明）。
%
% OUTPUT:  a : solution including boundary values
%          Q : reaction force vector
%              dim(a)=dim(Q)= nd x 1, nd : number of dof's
%-------------------------------------------------------------
% 缓存使用说明（重要，请先阅读）：
%  是否复用上次分解完全由第 4 个参数 reuse 决定（缺省 false）：
%     solveq(K,f) 或 solveq(K,f,bc)       → 重新分解并更新缓存（安全）
%     solveq(K,f,bc,true)                 → 复用缓存（O(1)，仅三角回代）
%
%  【启用缓存的条件（何时可传 reuse=true）】：
%    本次 K 与 bc 必须和最近一次 solveq 调用完全一致。
%    典型用法——主问题+伴随问题共用同一 K、同一 BC：
%        U   = solveq(K,F,BC);           % 主问题：分解并缓存
%        lam = solveq(K,psd,BC,true);    % 伴随问题：复用（O(1)）
%    【何时不能传 reuse=true】：K 或 bc 一旦改变（进入新迭代、换网格/边界），
%    必须先不带 reuse 调用一次重新分解，之后才能再复用。
%    缓存保留最近一次 {Cholesky 因子 L, 自由自由度集合, 置换 p}；`clear solveq` 可清空释放内存。
%
%  内存要点（3D 网格）：列主序编号下刚度矩阵 z 方向带宽极大，直接 chol 会产生灾难性
%  填充（如 60x60x20 → ~37GB）；因此分解前先用 dissect 做嵌套剖分填充最小化排序，
%  再对 Kr(p,p) 做稀疏 Cholesky，置换 p 一并缓存。
%
%  BC 的两种情形（使用不同式子）：
%    - 空矩阵 或 第 2 列全为 0（齐次边界）：直接解 d(free) = L'\(L\f(free))，d(pdof)=0；
%    - 第 2 列含非零给定位移（位移边界）：右端加已知位移项，
%      d(free) = L'\(L\(f(free) - K(free,pdof)*dp))，d(pdof)=dp。
%
%  Q (reaction forces) 仅在请求第 2 个输出时计算。
%  注意：K 必须对称正定（FE 刚度矩阵天然满足），chol 对其有效。
%  依赖：dissect（R2017a+，Graph and Network Algorithms 工具箱）。
%-------------------------------------------------------------

% LAST MODIFIED: H Yingvar   2026-08-13
% Changes:
% 1. Optimize memory usage and computation speed when solving very large matrices. Principle: use sparse `Cholesky` decomposition + `dissect` nested dissection ordering.
% 2. Add a reuse parameter; when [K] and [BC] are the same, allow reuse of cached data.
%-------------------------------------------------------------

  if nargin < 4, reuse = false; end                     % 缺省：不启用缓存
  if nargin < 3 || isempty(bc), bc = zeros(0, 2); end   % 无/空 bc → 齐次边界

  persistent Lcache freecache pcache                    % 缓存：Cholesky 因子 L、自由自由度集合、排序置换 p

  pdof  = bc(:, 1);
  dp    = bc(:, 2);
  homog = isempty(pdof) || all(dp == 0);                % 空或零 BC = 齐次边界
  nd    = size(K, 1);

  if reuse && ~isempty(Lcache)
      % ---- 用户指定复用：直接用缓存的 L/自由集合/置换（跳过 setdiff、排序与重分解）----
      L    = Lcache;
      free = freecache;
      p    = pcache;
  else
      % ---- 重新分解，并更新缓存 ----
      free = setdiff((1:nd)', pdof);
      Kr   = K(free, free);                             % 缩减刚度矩阵
      p    = dissect(Kr);                               % 嵌套剖分填充最小化排序（3D 网格必需，否则填充灾难）
      L    = chol(Kr(p, p), 'lower');                   % 排序后的稀疏 Cholesky 分解
      Lcache    = L;
      freecache = free;
      pcache    = p;
  end

  d = zeros(nd, 1);
  if homog
      % ---- 齐次边界（空或零 BC）：无需给定位移项 ----
      rhs = f(free);
  else
      % ---- 给定位移边界（非零 dp）：把已知位移移到右端 ----
      rhs = f(free) - K(free, pdof) * dp;
      d(pdof) = dp;
  end
  % 在置换空间求解：Kr(p,p)*x(p) = rhs(p) → x(p) = L'\(L\rhs(p))，再按 p 散回
  sol      = zeros(numel(free), 1);
  sol(p)   = L' \ (L \ rhs(p));
  d(free)  = sol;

  if nargout >= 2
      Q = K * d - f;                                     % reaction forces (only if asked)
  end

%--------------------------end--------------------------------
