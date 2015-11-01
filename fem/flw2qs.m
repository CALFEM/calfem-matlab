function [es,et]=flw2qs(ex,ey,ep,D,ed,eq)
% [es,et]=flw2qs(ex,ey,ep,D,ed,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute flows or corresponding quantities in the
%  quadrilateral field element.
% 
% INPUT:  ex = [x1 x2 x3 x4]
%         ey = [y1 y2 y3 y4]         element coordinates
%                             
%         ep = [t]                   element thickness 
%
%         D =  [kxx kxy;
%               kyx kyy]             constitutive matrix
%
%         ed =[u1 u2 u3 u4
%              .. .. .. ..]          u1,u2,u3,u4: nodal values
% 
%         eq                         heat supply per unit volume
%
% OUTPUT: es=[qx qy ] 
%             ... ..]                element flows
%
%         et=[gx gy ]
%             ... ..]                element gradients
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  K=zeros(5);  f=zeros(5,1);

  xm=sum(ex)/4; ym=sum(ey)/4;

  if nargin==5, q=0; else, q=eq; end

  En=[1 1 2 5; 2 2 3 5; 3 3 4 5; 4 4 1 5];
  ex1=[ex(1) ex(2) xm];
  ey1=[ey(1) ey(2) ym];
  ex2=[ex(2) ex(3) xm];
  ey2=[ey(2) ey(3) ym];
  ex3=[ex(3) ex(4) xm];
  ey3=[ey(3) ey(4) ym];
  ex4=[ex(4) ex(1) xm];
  ey4=[ey(4) ey(1) ym];

  [k1,f1]=flw2te(ex1,ey1,ep,D,q);
  [K,f]=assem(En(1,:),K,k1,f,f1);
  [k1,f1]=flw2te(ex2,ey2,ep,D,q);
  [K,f]=assem(En(2,:),K,k1,f,f1);
  [k1,f1]=flw2te(ex3,ey3,ep,D,q);
  [K,f]=assem(En(3,:),K,k1,f,f1);
  [k1,f1]=flw2te(ex4,ey4,ep,D,q);
  [K,f]=assem(En(4,:),K,k1,f,f1);

  [ni nj]=size(ed);

  a(5,ni)=0.;
  for i=1:ni
    a(:,i)=solveq(K,f,[[1:4]',ed(i,1:4)']);
  end

  [s1,t1]=flw2ts(ex1,ey1,D,[a(En(1,2:4),:)']);
  [s2,t2]=flw2ts(ex2,ey2,D,[a(En(2,2:4),:)']);
  [s3,t3]=flw2ts(ex3,ey3,D,[a(En(3,2:4),:)']);
  [s4,t4]=flw2ts(ex4,ey4,D,[a(En(4,2:4),:)']);

  es=(s1+s2+s3+s4)/4;
  et=(t1+t2+t3+t4)/4;
%--------------------------end--------------------------------
