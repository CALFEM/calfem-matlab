function [es,et,eci]=soli8s(ex,ey,ez,ep,D,ed)
% [es,et,eci]=soli8s2(ex,ey,ez,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for an
%  8 node (brick) isoparametric element.
%  ep=[1] (ir=1) evaluates at the element CENTRE point.
%
% INPUT:  ex = [x1 x2 x3 ... x8]   1 x 8 node x-coordinates
%         ey = [y1 y2 y3 ... y8]   1 x 8 node y-coordinates
%         ez = [z1 z2 z3 ... z8]   1 x 8 node z-coordinates
%         Each row corresponds to one element. If ex/ey/ez have a
%         single row, all elements share the same geometry (the
%         computation is then fully vectorized, with no loop).
%
%         ep = [Ir]                Ir: integration rule
%                                  Ir = 1 -> 1x1x1  (centre point)
%                                  Ir = 2 -> 2x2x2  (8 gauss points)
%                                  Ir = 3 -> 3x3x3  (27 gauss points)
%
%         D                        constitutive matrix (6 x 6)
%
%         ed = [u1 v1 w1 ... u8 v8 w8;   element displacement
%               ....................]    one row for each element
%              Node DOFs are interleaved: [u1 v1 w1 u2 v2 w2 ...].
%
% OUTPUT: es = [ sigx sigy sigz sigxy sigyz sigxz ;   stress, one row
%               ...............................   ]   per integration
%              point (element-major, point-minor order).
%         et = [ epsx epsy epsz gxy gyz gxz ; ... ]  strain, same row
%              ordering as es; shear entries are engineering strains.
%         eci= [ x y z ; ... ]  integration point coordinates.
%              Single element: ngp x 3. Multiple elements follow the
%              CALFEM convention of horizontal concatenation, giving
%              ngp x (3*rowed) (identical to the standard soli8s).
%-------------------------------------------------------------

% LAST MODIFIED: H Yingvar   2026-08-10
% Change Log:
% 1. Add case `ir = 1` for 1x1x1 integration point (centre point) with weight 8.
% 2. Adjust the cut components to the correct order： [sigxy sigyz sigxz].
% 3. Removed the integration-point for loop and compute all integration points in one vectorized pass (any ngp). 
% 4. All 3x3 Jacobian inverses are obtained at once via the adjugate (cofactor) formula. 
% 5. When all elements share the same geometry (rowex==1) the function is fully loop-free; otherwise (rowex>1) only the element loop is kept.

%-------------------------------------------------------------
  ir=ep(1);  ngp=ir*ir*ir;   % ir: 每个方向的积分点数, ngp=ir^3 为积分点总数

%--------- gauss points --------------------------------------
% 构造三维高斯积分点坐标 gp 与权重 w (三列分别对应 xi,eta,zeta 方向)
if ir==1
  gp = [0 0 0];   % 1x1x1 (centre point)
  w = [2 2 2];    % weight 2*2*2=8
elseif ir==2     % 2x2x2 = 8 个积分点 (8节点六面体的标准完全积分)
  g1=0.577350269189626; w1=1;   % 2点高斯坐标 +-1/sqrt(3), 权重均为1
  gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
  gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
  gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
elseif ir==3     % 3x3x3 = 27 个积分点
  % 一维 3 点高斯: 坐标 -g1, 0, +g1; 对应权重 w1(5/9), w2(8/9), w1(5/9)
  % 三维 27 个点 = 三个方向坐标的张量积; 行序: xi 变化最快 -> eta -> zeta 最慢
  g1=0.774596669241483;            % 3点高斯坐标幅值: -0.775, 0, +0.775
  w1=0.555555555555555; w2=0.888888888888888;   % 3点高斯权重: 5/9 与 8/9

  coords=[-g1 0 g1];               % 一维坐标向量
  w1d   =[w1 w2 w1];               % 一维权重向量 (与 coords 一一对应)
  [xi,eta,zeta]    = ndgrid(coords);  % 3x3x3 张量网格
  [wxi,weta,wzeta] = ndgrid(w1d);
  gp = [xi(:) eta(:) zeta(:)];     % 27 x 3 积分点坐标
  w  = [wxi(:) weta(:) wzeta(:)];  % 27 x 3 各方向的一维权重
else
  disp('Used number of integration points not implemented');
  return
end

wp=w(:,1).*w(:,2).*w(:,3);
xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;  % 自然坐标分量; r2=3*ngp

%--------- shape functions -----------------------------------
% 8节点三线性形函数 N (ngp x 8): 每行对应一个积分点, 每列对应一个节点
N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

% dNr (3*ngp x 8): 形函数对自然坐标的导数, 每个积分点占 3 行
%   行 3i-2 = d/dxi, 行 3i-1 = d/deta, 行 3i = d/dzeta
dNr(1:3:r2,1)=-(1-eta).*(1-zet);    dNr(1:3:r2,2)= (1-eta).*(1-zet);
dNr(1:3:r2,3)= (1+eta).*(1-zet);    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
dNr(1:3:r2,5)=-(1-eta).*(1+zet);    dNr(1:3:r2,6)= (1-eta).*(1+zet);
dNr(1:3:r2,7)= (1+eta).*(1+zet);    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
dNr=dNr/8.;   % 统一除以 8

%--------- three dimensional case ----------------------------
rowed=size(ed,1);   % 位移向量 ed 的行数 = 待处理的单元个数
rowex=size(ex,1);   % 坐标矩阵的行数

  % 若 ex 只有 1 行, 表示所有单元共用同一套几何坐标 (ie 不递增); 否则每行一个单元
  if rowex==1 
    incie=0; 
  else 
    incie=1; 
  end

  % 预分配输出: es/et 每单元 ngp 个积分点 (row*ngp 行);
  % eci 兼容旧版 [eci X] 横向拼接约定: ngp 行 x (3*rowed) 列
  es = zeros(rowed*ngp,6);
  et = zeros(rowed*ngp,6);
  eci= zeros(ngp, 3*rowed);

  % 将 dNr 重排为 3 x 8 x ngp: 第 k 页 dNr3(:,:,k) = 积分点 k 的 3x8 导数块
  dNr3 = permute(reshape(dNr,3,ngp,8),[1 3 2]);

  if rowex==1
    % ============ 所有单元共用同一几何: 完全向量化, 无任何 for 循环 ============
    co = [ex; ey; ez];                     % 3 x 8 单元节点坐标
    eci_all = N*co';                       % 积分点插值坐标 (ngp x 3), 各单元相同

    JT = dNr*co';                          % 3*ngp x 3
    [d1,d2,d3] = jacDerivs(JT, dNr3, ngp); % dN/dx, dN/dy, dN/dz (各 8 x ngp)

    U = ed(:,1:3:22);  V = ed(:,2:3:23);  W = ed(:,3:3:24);  % rowed x 8 三向位移
    % 六个应变分量 (rowed x ngp), 顺序 [epsx epsy epsz gxy gyz gxz]
    ee6 = cat(3, U*d1, V*d2, W*d3, U*d2+V*d1, V*d3+W*d2, U*d3+W*d1);
    % 重排为 (rowed*ngp) x 6, 行序按 单元->积分点 (与逐点循环顺序一致)
    et = reshape(permute(ee6,[2 1 3]), rowed*ngp, 6);
    es = et*D';                            % 应力 = 应变 * D'
    eci= repmat(eci_all, 1, rowed);        % ngp x (3*rowed), 与旧版横向拼接一致
  else
    % ============ 各单元几何不同: 只保留单元循环, 积分点已全部向量化 ============
    ie=1;                                  % 当前单元坐标索引
    for ied=1:rowed
      co = [ex(ie,:); ey(ie,:); ez(ie,:)];                 % 3 x 8
      eci(:,(ied-1)*3+1:ied*3) = N*co';                    % 横向拼接 (ngp x 3) 块

      JT = dNr*co';                                        % 3*ngp x 3
      [d1,d2,d3] = jacDerivs(JT, dNr3, ngp);               % dN/dx, dN/dy, dN/dz (8 x ngp)

      u = ed(ied,1:3:22)';  v = ed(ied,2:3:23)';  w = ed(ied,3:3:24)';  % 8 节点三向位移
      % 应变 (6 x ngp): 与原先 B 矩阵 [epsx epsy epsz gxy gyz gxz] 完全一致
      ee = [ u'*d1; v'*d2; w'*d3; u'*d2+v'*d1; v'*d3+w'*d2; u'*d3+w'*d1 ];
      et((ied-1)*ngp+1:ied*ngp,:) = ee';                   % 每积分点一行
      es((ied-1)*ngp+1:ied*ngp,:) = (D*ee)';               % 应力 = 本构矩阵 D * 应变

      ie=ie+incie;    % 移动到下一个单元的坐标索引
    end
  end
%--------------------------end--------------------------------

% ==================== 局部函数: 向量化求形函数物理导数 ====================
function [d1,d2,d3] = jacDerivs(JT, dNr3, ngp)
% 一次性求出所有积分点的形函数对物理坐标的导数 (dN/dx, dN/dy, dN/dz)
%   JT  : 3*ngp x 3, 每连续 3 行构成一个积分点的 3x3 雅可比矩阵
%   dNr3: 3 x 8 x ngp, 形函数对自然坐标的导数 (第 k 页为积分点 k)
%   d1,d2,d3: 8 x ngp, 分别为 dN/dx, dN/dy, dN/dz
  JT3 = permute(reshape(JT,3,ngp,3),[1 3 2]);   % 3 x 3 x ngp: 各积分点的雅可比矩阵

  % 3x3 行列式的余子式 (1 x 1 x ngp), 供伴随矩阵求逆使用
  c11=JT3(2,2,:).*JT3(3,3,:)-JT3(2,3,:).*JT3(3,2,:);
  c12=-(JT3(2,1,:).*JT3(3,3,:)-JT3(2,3,:).*JT3(3,1,:));
  c13= JT3(2,1,:).*JT3(3,2,:)-JT3(2,2,:).*JT3(3,1,:);
  c21=-(JT3(1,2,:).*JT3(3,3,:)-JT3(1,3,:).*JT3(3,2,:));
  c22= JT3(1,1,:).*JT3(3,3,:)-JT3(1,3,:).*JT3(3,1,:);
  c23=-(JT3(1,1,:).*JT3(3,2,:)-JT3(1,2,:).*JT3(3,1,:));
  c31= JT3(1,2,:).*JT3(2,3,:)-JT3(1,3,:).*JT3(2,2,:);
  c32=-(JT3(1,1,:).*JT3(2,3,:)-JT3(1,3,:).*JT3(2,1,:));
  c33= JT3(1,1,:).*JT3(2,2,:)-JT3(1,2,:).*JT3(2,1,:);

  detJ3 = JT3(1,1,:).*c11 + JT3(1,2,:).*c12 + JT3(1,3,:).*c13;  % 1 x 1 x ngp
  if any(detJ3(:) < 10*eps)
    disp('Jacobideterminant equal or less than zero!')
  end

  % dNx = inv(J)*dNr = adj(J)*dNr/det, 逐行展开, 所有积分点一次算完
  dNx3 = zeros(3,8,ngp);
  dNx3(1,:,:) = (c11.*dNr3(1,:,:) + c21.*dNr3(2,:,:) + c31.*dNr3(3,:,:))./detJ3;
  dNx3(2,:,:) = (c12.*dNr3(1,:,:) + c22.*dNr3(2,:,:) + c32.*dNr3(3,:,:))./detJ3;
  dNx3(3,:,:) = (c13.*dNr3(1,:,:) + c23.*dNr3(2,:,:) + c33.*dNr3(3,:,:))./detJ3;

  % 注意: 不能直接用 squeeze —— ngp=1 时 1x8x1 会被 squeeze 成 1x8 行向量,
  % 必须用 reshape 显式保持 8 x ngp 的列方向
  d1 = reshape(dNx3(1,:,:), 8, ngp);   % 8 x ngp: dN/dx
  d2 = reshape(dNx3(2,:,:), 8, ngp);   % 8 x ngp: dN/dy
  d3 = reshape(dNx3(3,:,:), 8, ngp);   % 8 x ngp: dN/dz
