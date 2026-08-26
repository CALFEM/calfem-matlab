%% test_soli8s2.m
% soli8s2 —— 8 节点六面体单元应力/应变恢复函数的单元测试
%
% 运行方式:
%   1) 在 fem 目录下的命令行输入:   test_soli8s2
%   2) 或打开本文件后点编辑器"运行"(F5)
%
% 结果判读:
%   - 每一项打印相对误差和 PASS/FAIL;
%   - 常应变场(仿射位移场)对三线性单元应当"精确"重现,
%     相对误差应在 ~1e-15(机器精度)量级;
%   - 末尾输出 "ALL TESTS PASSED" 即全部通过; 若出现 FAIL 或
%     误差明显大于 1e-9, 说明实现有误。
%   - 测试覆盖: 单单元 ir=1/2/3、共享几何多单元、几何不同多单元、
%     应力-应变一致性、ir=1 中心点坐标。

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
nPass = 0; nFail = 0;
relerr = @(a,b) max(max(abs(a-b))) / max(1, max(max(abs(a))));

%% =====================================================================
%  Test 1: 单单元, ir=1/2/3 —— 常应变场精确重现
%  常应变场 u=eps_x*x+(gxy/2)*y+(gxz/2)*z, ... (仿射) 三线性单元精确重现
%% =====================================================================
fprintf('==== Test 1: 单单元常应变场 (ir=1,2,3) ====\n');
strain0 = [0.001; -0.0005; 0.0008; 0.0006; -0.0004; 0.0003];  % [epsx epsy epsz gxy gyz gxz]
for ir = 1:3
    ed = zeros(1,24);                                    % [u1 v1 w1 ... u8 v8 w8]
    ed(1:3:22) = strain0(1)*x + (strain0(4)/2)*y + (strain0(6)/2)*z;   % u
    ed(2:3:23) = (strain0(4)/2)*x + strain0(2)*y + (strain0(5)/2)*z;   % v
    ed(3:3:24) = (strain0(6)/2)*x + (strain0(5)/2)*y + strain0(3)*z;   % w

    [es,et,eci] = soli8s2(x,y,z,[ir],D,ed);
    ngp = ir^3;
    e1 = relerr(et, repmat(strain0', ngp, 1));            % 应变 vs 期望
    e2 = relerr(es, repmat((D*strain0)', ngp, 1));        % 应力 vs 期望
    ok = (e1<tol && e2<tol && isequal(size(es),[ngp 6]) && isequal(size(eci),[ngp 3]));
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  应变相对差=%.2e  应力相对差=%.2e  尺寸es=%s eci=%s  -> %s\n', ...
        ir, e1, e2, mat2str(size(es)), mat2str(size(eci)), status(ok));
end

%% =====================================================================
%  Test 2: ir=1 中心点 —— 积分点物理坐标应为单元"形心"(节点坐标平均)
%% =====================================================================
fprintf('==== Test 2: ir=1 中心点坐标 ====\n');
ed1 = zeros(1,24); ed1(1:3:22) = 0.001*x;                % 任意单向场
[~,~,eci1] = soli8s2(x,y,z,[1],D,ed1);
center = mean([x;y;z],2)';                               % 期望中心点 [0 0 0]
ok = all(abs(eci1 - center) < tol);
if ok, nPass=nPass+1; else, nFail=nFail+1; end
fprintf('  eci=[%.6g %.6g %.6g]  期望=[%.6g %.6g %.6g]  -> %s\n', ...
    eci1, center, status(ok));

%% =====================================================================
%  Test 3: 共享几何多单元 (rowex==1, rowed>1) —— 每单元应变放大 f 倍
%% =====================================================================
fprintf('==== Test 3: 共享几何多单元 (rowex==1) ====\n');
N = 5;
edm = zeros(N,24);
for k = 1:N
    f = k;
    edm(k,1:3:22) = f*(strain0(1)*x + (strain0(4)/2)*y + (strain0(6)/2)*z);
    edm(k,2:3:23) = f*((strain0(4)/2)*x + strain0(2)*y + (strain0(5)/2)*z);
    edm(k,3:3:24) = f*((strain0(6)/2)*x + (strain0(5)/2)*y + strain0(3)*z);
end
for ir = 1:3
    ngp = ir^3;
    [es,et,eci] = soli8s2(x,y,z,[ir],D,edm);
    % 期望: 行序 单元->积分点, 单元 k 的应变 = k*strain0
    et_exp = repelem((1:N)', ngp, 1) .* repmat(strain0', N*ngp, 1);
    es_exp = et_exp*D';
    e1 = relerr(et, et_exp); e2 = relerr(es, es_exp);
    sz = isequal(size(es),[N*ngp 6]) && isequal(size(et),[N*ngp 6]) ...
         && isequal(size(eci),[ngp 3*N]);               % eci 保持 CALFEM 横向拼接形状
    ok = (e1<tol && e2<tol && sz);
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  应变相对差=%.2e  应力相对差=%.2e  尺寸es=%s eci=%s  -> %s\n', ...
        ir, e1, e2, mat2str(size(es)), mat2str(size(eci)), status(ok));
end

%% =====================================================================
%  Test 4: 几何不同多单元 (rowex>1) —— 缩放/平移的立方体, 仍应精确重现常应变
%% =====================================================================
fprintf('==== Test 4: 几何不同多单元 (rowex>1) ====\n');
N2 = 4;
exm = zeros(N2,8); eym = zeros(N2,8); ezm = zeros(N2,8);
for k = 1:N2
    s  = 0.8 + 0.3*k;                                  % 缩放
    sh = 0.2*k;                                        % 平移
    exm(k,:) = s*x + sh; eym(k,:) = s*y + sh; ezm(k,:) = s*z + sh;
end
% 各单元构造相同的常应变场: 位移必须用"该单元的实际物理坐标"
edm2 = zeros(N2,24);
for k = 1:N2
    edm2(k,1:3:22) = strain0(1)*exm(k,:) + (strain0(4)/2)*eym(k,:) + (strain0(6)/2)*ezm(k,:);
    edm2(k,2:3:23) = (strain0(4)/2)*exm(k,:) + strain0(2)*eym(k,:) + (strain0(5)/2)*ezm(k,:);
    edm2(k,3:3:24) = (strain0(6)/2)*exm(k,:) + (strain0(5)/2)*eym(k,:) + strain0(3)*ezm(k,:);
end
for ir = 1:3
    ngp = ir^3;
    [es,et,eci] = soli8s2(exm,eym,ezm,[ir],D,edm2);
    et_exp = repmat(strain0', N2*ngp, 1);
    es_exp = et_exp*D';
    e1 = relerr(et, et_exp); e2 = relerr(es, es_exp);
    sz = isequal(size(es),[N2*ngp 6]) && isequal(size(eci),[ngp 3*N2]);
    ok = (e1<tol && e2<tol && sz);
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  应变相对差=%.2e  应力相对差=%.2e  尺寸es=%s eci=%s  -> %s\n', ...
        ir, e1, e2, mat2str(size(es)), mat2str(size(eci)), status(ok));
end

%% =====================================================================
%  Test 5: 一致性 —— 任意位移场下都应满足 es = et*D' (顺序/公式不串位)
%% =====================================================================
fprintf('==== Test 5: 应力-应变一致性 es = et*D'' ====\n');
rng(1);
edr = rand(1,24)*1e-3;
for ir = 1:3
    [es,et] = soli8s2(x,y,z,[ir],D,edr);
    e = relerr(es, et*D');
    ok = e<tol;
    if ok, nPass=nPass+1; else, nFail=nFail+1; end
    fprintf('  ir=%d  相对差=%.2e  -> %s\n', ir, e, status(ok));
end

%% =====================================================================
fprintf('\n================================================\n');
if nFail==0
    fprintf('ALL TESTS PASSED  (%d/%d)\n', nPass, nPass+nFail);
else
    fprintf('FAILED: %d/%d 项未通过, 请检查 soli8s2.m\n', nFail, nPass+nFail);
end

function s = status(ok)
% 简单 PASS/FAIL 标记
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
