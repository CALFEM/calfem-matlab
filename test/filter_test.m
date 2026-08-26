nelx = 10;
nely = 10;
rmin = 2;
% L形设计域
passive = zeros(nely, nelx);
passive(1:round(nely*0.6), round(nelx*0.4)+1:end) = 1; % 挖去部分

% 创建圆形设计域（中心在网格中心）
% [x, y] = meshgrid(1:nelx, 1:nely);
% center_x = nelx/2; center_y = nely/2;
% radius = min(nelx, nely)/2.5;  % 半径
% passive = sqrt((x-center_x).^2 + (y-center_y).^2) > radius;

% % 内外半径定义环形
% inner_radius = min(nelx, nely)/5;
% outer_radius = min(nelx, nely)/2.5;
% passive = ~(sqrt((x-center_x).^2 + (y-center_y).^2) <= outer_radius & ...
%            sqrt((x-center_x).^2 + (y-center_y).^2) >= inner_radius);
% 
% % 等边三角形（顶点在网格顶部和两侧）
% passive = true(nely, nelx);
% for i = 1:nelx
%     height = nely*(1 - abs(i - nelx/2)/(nelx/2));
%     passive(1:round(height), i) = false;  % 三角形内部为设计域
% end


[Hs, Hker] = elemental_filter_prepare(nelx, nely, rmin, passive);
Hs(logical(passive)) = 0;
% Hs = full(reshape(Hs,nelx,nely));
imagesc(Hs);colorbar;axis image;
% Hs_full = full(reshape(Hs,nelx,nely));
% [dy,dx] = meshgrid(-ceil(rmin)+1: ceil(rmin)-1,-ceil(rmin)+1: ceil(rmin)-1);
% h = max(0,rmin-sqrt(dx.^2+dy.^2));
% alpha_i = rand(nelx,nely);
% tic;
% alpha_i1 = H1 * alpha_i(:) ./ Hs(:);
% t1 = toc;
% disp(t1);
% tic;
% alpha_i2 = conv2(alpha_i,h,"same") ./ Hs;
% t2 = toc;
% disp(t2);
% disp(["conv2 比传统快:",(t1-t2)/t1 * 100, "%"]);
% % 得出结论，conv2更快
% Hs2 = conv2(~passive,h,"same");
% Hs2(passive == 1) = 0;
% figure;
% imagesc(Hs2);colorbar; axis image;

% [wsum,ker, w1] = nodal_filter_prepare(nelx, nely, rmin, passive);
% wsum(logical(passive)) = 0;
% imagesc(wsum);colorbar; axis image;