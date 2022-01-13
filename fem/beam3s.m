function [es,edi,eci]=beam3s(ex,ey,ez,eo,ep,ed,eq,n)
%  es=beam3s(ex,ey,ez,eo,ep,ed)
%  es=beam3s(ex,ey,ez,eo,ep,ed,eq)
% [es,edi]=beam3s(ex,ey,ez,eo,ep,ed,eq,n)
% [es,edi,eci]=beam3s(ex,ey,ez,eo,ep,ed,eq,n)
%-----------------------------------------------------------------------
%    PURPOSE:
%     Calculate the variation of the section forces and displacements
%     along a three-dimensional beam element
%
%    INPUT:  ex = [x1 x2]       
%            ey = [y1 y2]      
%            ez = [z1 z2]             node coordinates
%
%            eo = [xz yz zz]          orientation of local z-axis  
%
%
%            ep = [E G A Iy Iz Kv]    element properties:
%                                     E = Young's modulus
%                                     G = Shear modulus 
%                                     A = the cross section area
%                                     Iy= the moment of inertia, local y-axis
%                                     Iz= the moment of inertia, local z-axis
%                                     Kv= Saint-Venant's torsion constant
%
%            ed =                     the element displacement vector from the
%                                     global coordinate system
%
%            eq = [qx qy qz qw]       the distributed axial, transversal and
%                                     torsional loads
%
%             n =                     the number of points in which displacements
%                                     and section forces are to be computed
%
%   OUTPUT:
%
%    es = [N1 Vy1 Vz1 T1 My1 Mz1;      section forces in n points
%          N2 Vy2 Vz2 T2 My2 Mz2;      along the local x-axis
%          .. ... ... .. ... ...;
%          Nn Vyn Vzn Tn Myn Mzn]
%
%    edi = [u1 v1 w1 fi1;              displacements in n points 
%           u2 v2 w2 fi2;              along the local x-axis
%           .. .. .. ...;
%           un vn wn fin]
%
%    eci = [x1;                        evaluation points on the 
%           x2;                        local x-axis 
%            .
%           xn]
%
%-----------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2021-09-06
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-----------------------------------------------------------------------
  E=ep(1); Gs=ep(2);
  A=ep(3); Iy=ep(4); Iz=ep(5); Kv=ep(6);
  DEA=E*A; DEIz=E*Iz; DEIy=E*Iy; DGK=Gs*Kv;

  qX=0; qY=0; qZ=0; qW=0; 
  if nargin>6; qX=eq(1); qY=eq(2); qZ=eq(3); qW=eq(4); end

  ne=2;        if nargin>7;  ne=n; end;
  
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  dz=ez(2)-ez(1);
  L=sqrt(dx*dx+dy*dy+dz*dz);
  n1=[dx dy dz]/L;
  lc=sqrt(eo*eo'); 
  n3=eo/lc;
   
  n2(1)= n3(2)*n1(3)-n3(3)*n1(2);
  n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
  n2(3)= n3(1)*n1(2)-n1(1)*n3(2);
%
  An=[n1; n2; n3];
%
  G=[  An     zeros(3) zeros(3) zeros(3);
     zeros(3)   An     zeros(3) zeros(3);
     zeros(3) zeros(3)   An     zeros(3);
     zeros(3) zeros(3) zeros(3)   An    ];

  edl=G*ed';

  a1=[edl(1); edl(7)];
  C1=[ 1   0;
      -1/L 1/L];
  C1a=C1*a1;
  
  a2=[edl(2); edl(6); edl(8); edl(12)];
  C2=[1       0    0       0;
      0       1    0       0;
     -3/(L^2) -2/L 3/(L^2) -1/L;
     2/(L^3) 1/(L^2) -2/(L^3) 1/(L^2)];
  C2a=C2*a2;

  a3=[edl(3); -edl(5); edl(9); -edl(11)];
  C3=C2;
  C3a=C3*a3;
  
  a4=[edl(4); edl(10)];
  C4=C1;
  C4a=C4*a4;
  
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));

  u=[one X]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u-(X.^2-L*X)*qX/(2*DEA);
    du=du-(2*X-L)*qX/(2*DEA);
  end; 
  v=[one X X.^2 X.^3]*C2a;
  d2v=[zero zero 2*one 6*X]*C2a;
  d3v=[zero zero zero 6*one]*C2a;
  if DEIz~=0;
    v=v+(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEIz);
    d2v=d2v+(6*X.^2-6*L*X+L^2)*qY/(12*DEIz);
    d3v=d3v+(2*X-L)*qY/(2*DEIz);
  end;
  w=[one X X.^2 X.^3]*C3a;
  d2w=[zero zero 2*one 6*X]*C3a;
  d3w=[zero zero zero 6*one]*C3a;
  if DEIy~=0;
    w=w+(X.^4-2*L*X.^3+L^2*X.^2)*qZ/(24*DEIy);
    d2w=d2w+(6*X.^2-6*L*X+L^2)*qZ/(12*DEIy);
    d3w=d3w +(2*X-L)*qZ/(2*DEIy);
  end;
  fi=[one X]*C4a;
  dfi=[zero one]*C4a;
  if DGK~=0;
    fi=fi-(X.^2-L*X)*qW/(2*DGK);
    dfi=dfi-(2*X-L)*qW/(2*DGK);
  end;  
  
  N=DEA*du; Mz=DEIz*d2v; Vy=-DEIz*d3v;  My=-DEIy*d2w; Vz=-DEIy*d3w; T=DGK*dfi; 
  es=[N Vy Vz T My Mz];
  edi=[u v w fi];
  eci=X;
%----------------------------------end----------------------------------
