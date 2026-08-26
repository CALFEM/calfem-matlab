%% test_soli8s_fields.m
% soli8s —— 用"不同位移场"补充验证 (与 test_soli8s.m 互补)
%
% 运行方式:
%   在 fem 目录命令行输入:  test_soli8s_fields
%   或打开本文件后点"运行"(F5)
%
% 本文件用更丰富的位移场验证结果是否正确:
%   F0  常应变(仿射)场 —— 三线性单元应精确重现 (分片检验)
%   F1  一般三线性场 (应变随位置变化, 逐点"精确"重现)
%       —— 可抓住"常应变场验不出来"的 bug, 如积分点应变错位/乱序;
%   F2  刚体平移 + 刚体转动 (期望应变/应力 = 0)
%       —— 检验位移梯度的反对称(旋转)部分不会漏进应变;
%   F3  二次场 u = c*x^2 (不在三线性单元空间内): 网格加密后应变
%       应收敛到解析梯度 epsx = 2cx —— 端到端精度/收敛性验证。
%
% 判读: 每项打印相对误差(或绝对值)与 PASS/FAIL, 末尾汇总;
%       F0/F1 误差应在 ~1e-16(机器精度), F2 应 ~0, F3 误差应随网格
%       加密(n=1,2,4)单调减小。

clear; clc;

%% ---------- 材料(钢, SI)与 6x6 本构矩阵 D ----------
E = 210e9;  nu = 0.3;
c  = E/((1+nu)*(1-2*nu));
D  = c*[1-nu nu nu 0 0 0;
        nu 1-nu nu 0 0 0;
        nu nu 1-nu 0 0 0;
        0 0 0 (1-2*nu)/2 0 0;
        0 0 0 0 (1-2*nu)/2 0;
        0 0 0 0 0 (1-2*nu)/2];

%% ---------- 标准 8 节点立方体单元节点坐标 (-1..1) ----------
x = [-1  1  1 -1 -1  1  1 -1];
y = [-1 -1  1  1 -1 -1  1  1];
z = [-1 -1 -1 -1  1  1  1  1];

tol   = 1e-9;               % 相对误差容差
tol0  = 1e-12;              % F2 绝对容差(期望精确为 0)
nPass = 0; nFail = 0;
relerr = @(a,b) max(max(abs(a-b))) / max(1, max(max(abs(a))));

%% =====================================================================
%  F0: 常应变(仿射)场 —— 三线性单元应精确重现 (分片检验/patch test)
%  u=eps_x*x+(gxy/2)*y+(gxz/2)*z, ... 应变处处等于 strain0
%% =====================================================================
fprintf('==== F0: 常应变场 (ir=1,2,3) ====\n');
strain0 = [0.001; -0.0005; 0.0008; 0.0006; -0.0004; 0.0003];  % [epsx epsy epsz gxy gyz gxz]
for ir = 1:3
    ed = zeros(1,24);
    ed(1:3:22) = strain0(1)*x + (strain0(4)/2)*y + (strain0(6)/2)*z;   % u
    ed(2:3:23) = (strain0(4)/2)*x + strain0(2)*y + (strain0(5)/2)*z;   % v
    ed(3:3:24) = (strain0(6)/2)*x + (strain0(5)/2)*y + strain0(3)*z;   % w
    [es,et] = soli8s(x,y,z,[ir],D,ed);
    ngp = ir^3;
    e1 = relerr(et, repmat(strain0', ngp, 1));            % 应变 vs 期望
    e2 = relerr(es, repmat((D*strain0)', ngp, 1));        % 应力 vs 期望
    ok = (e1<tol && e2<tol);
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  应变相对差=%.2e  应力相对差=%.2e  -> %s\n', ...
        ir, e1, e2, status(ok));
end

%% =====================================================================
%  F1: 一般三线性场 —— 应变随位置变化, 每个积分点应"精确"重现解析值
%  场(均在单元空间内):  u=a1*x+a2*x*y+a3*x*z, v=b1*y+b2*x*y, w=g1*z+g2*y*z
%  解析应变在积分点 k 处(物理坐标取自 eci):
%   epsx=a1+a2*y+a3*z, epsy=b1+b2*x, epsz=g1+g2*y,
%   gxy=a2*x+b2*y, gyz=g2*z, gxz=a3*x
%% =====================================================================
fprintf('==== F1: 一般三线性场 (逐点应变不同) ====\n');
a1=1e-3; a2=2e-4; a3=3e-4; b1=1.5e-3; b2=2.5e-4; g1=8e-4; g2=1.2e-4;
ed = zeros(1,24);
ed(1:3:22) = a1*x + a2*x.*y + a3*x.*z;   % u
ed(2:3:23) = b1*y + b2*x.*y;             % v
ed(3:3:24) = g1*z + g2*y.*z;             % w
for ir = 1:3
    [es,et,eci] = soli8s(x,y,z,[ir],D,ed);
    et_exp = [a1 + a2*eci(:,2) + a3*eci(:,3), ...   % epsx
              b1 + b2*eci(:,1), ...                 % epsy
              g1 + g2*eci(:,2), ...                 % epsz
              a2*eci(:,1) + b2*eci(:,2), ...        % gxy
              g2*eci(:,3), ...                      % gyz
              a3*eci(:,1)];                         % gxz
    es_exp = et_exp*D';
    e1 = relerr(et, et_exp); e2 = relerr(es, es_exp);
    ok = (e1<tol && e2<tol);
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  应变相对差=%.2e  应力相对差=%.2e  -> %s\n', ...
        ir, e1, e2, status(ok));
end

%% =====================================================================
%  F2: 刚体平移 + 刚体转动 —— 应变/应力应为 0
%  平移: u=[a b c] 常数;  转动: u = omega × r
%% =====================================================================
fprintf('==== F2: 刚体位移/转动 -> 零应变 ====\n');
% 刚体平移
edT = repmat([1e-3 2e-3 -1.5e-3], 1, 8);
% 刚体转动 u = ω × r
wx=1e-4; wy=-2e-4; wz=3e-4;
edR = zeros(1,24);
edR(1:3:22) = -wz*y + wy*z;   % u = -θz*y + θy*z
edR(2:3:23) =  wz*x - wx*z;   % v =  θz*x - θx*z
edR(3:3:24) = -wy*x + wx*y;   % w = -θy*x + θx*y
for ir = 1:3
    for t = 1:2
        if t==1, ed = edT; name = '平移'; else, ed = edR; name = '转动'; end
        [es,et] = soli8s(x,y,z,[ir],D,ed);
        e_et = max(abs(et(:)));   % 应变: 应精确为 0(机器精度 ~1e-19)
        e_es = max(abs(es(:)));   % 应力 = D*应变, 底噪被 D(~1e11)放大成 ~1e-8 Pa
        ok = (e_et < 1e-12) && (e_es < 1e-6);   % 应力只检查"绝对底噪"量级
        if ok, nPass=nPass+1; else, nFail=nFail+1; end
        fprintf('  ir=%d %s  最大|应变|=%.2e  最大|应力|=%.2e Pa  -> %s\n', ...
            ir, name, e_et, e_es, status(ok));
    end
end

%% =====================================================================
%  F3: 二次场 u = c*x^2 (不在单元空间内) —— 网格收敛到解析梯度
%  解析应变 epsx = 2cx (随 x 线性变化), 其余为 0;
%  三线性单元不能精确重现, 但网格加密后应变应收敛到解析值。
%% =====================================================================
fprintf('==== F3: 二次场 u=c*x^2 -> 网格收敛 (ir=2 积分点) ====\n');
cq = 1e-3;
errN = zeros(1,4);
for n = [1 2 4]
    [exg,eyg,ezg,edg] = genQuadMesh(n, cq);
    ne = size(edg,1);
    et_all = zeros(ne*8,6); eci_all = zeros(ne*8,3);
    p = 0;
    for e = 1:ne
        [~, et1, eci1] = soli8s(exg(e,:), eyg(e,:), ezg(e,:), [2], D, edg(e,:));
        et_all(p+1:p+8,:) = et1;  eci_all(p+1:p+8,:) = eci1;  p = p+8;
    end
    exact = [2*cq*eci_all(:,1), zeros(size(eci_all,1),5)];   % epsx=2cx, 其余 0
    errN(n) = max(max(abs(et_all - exact)));
end
fprintf('  最大应变误差:  n=1: %.3e   n=2: %.3e   n=4: %.3e\n', ...
    errN(1), errN(2), errN(4));
ok = (errN(4) < errN(2) && errN(2) < errN(1));   % 网格加密误差应单调减小
if ok, nPass=nPass+1; else, nFail=nFail+1; end
fprintf('  收敛性(误差随加密单调减小) -> %s\n', status(ok));

%% =====================================================================
fprintf('\n================================================\n');
if nFail==0
    fprintf('ALL TESTS PASSED  (%d/%d)\n', nPass, nPass+nFail);
else
    fprintf('FAILED: %d/%d 项未通过, 请检查 soli8s.m\n', nFail, nPass+nFail);
end

function s = status(ok)
% 简单 PASS/FAIL 标记
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function [exg,eyg,ezg,edg] = genQuadMesh(n, cq)
% 生成 [-1,1]^3 上 n×n×n 的 8 节点六面体网格;
% 节点位移场 u = cq*x^2, v = w = 0 (二次场, 不在三线性单元空间内)
% 单元局部节点序与 soli8s 一致:
%   [(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)]
h  = 2/n;  Np = n+1;
X = zeros(Np,Np,Np); Y = X; Z = X;
for k = 0:n
    for j = 0:n
        for i = 0:n
            X(i+1,j+1,k+1) = -1 + h*i;
            Y(i+1,j+1,k+1) = -1 + h*j;
            Z(i+1,j+1,k+1) = -1 + h*k;
        end
    end
end
idx = @(i,j,k) i + Np*j + Np*Np*k + 1;    % 1-based 线性节点号
ne = n^3;
exg = zeros(ne,8); eyg = zeros(ne,8); ezg = zeros(ne,8); edg = zeros(ne,24);
e = 0;
for k = 0:n-1
    for j = 0:n-1
        for i = 0:n-1
            e = e+1;
            nlist = [idx(i,j,k),   idx(i+1,j,k),   idx(i+1,j+1,k), idx(i,j+1,k), ...
                     idx(i,j,k+1), idx(i+1,j,k+1), idx(i+1,j+1,k+1), idx(i,j+1,k+1)];
            exg(e,:) = X(nlist); eyg(e,:) = Y(nlist); ezg(e,:) = Z(nlist);
            edg(e,1:3:22) = cq*exg(e,:).^2;     % u = cq*x^2; v=w=0 保持 0
        end
    end
end
end
