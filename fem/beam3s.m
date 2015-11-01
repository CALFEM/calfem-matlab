function [es,edi,eci]=beam3s(ex,ey,ez,eo,ep,ed,eq,n)
%  es=beam3s(ex,ey,ez,eo,ep,ed)
%  es=beam3s(ex,ey,ez,eo,ep,ed,eq)
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
%    eci = [x1;
%           x2;
%            .
%            .
%           xn]              local x-coordinates of the evaluation
%                                      points
%-----------------------------------------------------------------------

% LAST MODIFIED: E Serrano     1995-09-21 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-----------------------------------------------------------------------
  if nargin<=7  n=2; end;
  if nargin>6
    qx=eq(1); qy=eq(2); qz=eq(3); qw=eq(4);
  else
    qx=0;qy=0;qz=0;qw=0;
  end 
%
  b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
  L=sqrt(b'*b);n1=b/L;
%
  lc=sqrt(eo*eo');n3=eo/lc;
%
  EIy=ep(1)*ep(4); EIz=ep(1)*ep(5);
  EA=ep(1)*ep(3); GKv=ep(2)*ep(6);
%   
  n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
  n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
  n2(3)=n3(1)*n1(2)-n1(1)*n3(2);
%
  An=[n1';
      n2;
      n3];
%
  G=[  An     zeros(3) zeros(3) zeros(3);
     zeros(3)   An     zeros(3) zeros(3);
     zeros(3) zeros(3)   An     zeros(3);
     zeros(3) zeros(3) zeros(3)   An    ];
%
  u=G*ed' - [0;                 % u is the local element dis-
             0;                 % placement vector minus the
             0;                 % particular solution to the  
             0;                 % beam's diff.eq:s
             0;
             0;
        -qx*L^2/2/EA;
         qy*L^4/24/EIz;
         qz*L^4/24/EIy;
        -qw*L^2/2/GKv;
        -qz*L^3/6/EIy;
         qy*L^3/6/EIz];
%
  C=[0 1   0    0  0 0    0     0   0 0 0 0;
     0 0   0    0  0 1    0     0   0 0 0 0;
     0 0   0    0  0 0    0     0   0 1 0 0;
     0 0   0    0  0 0    0     0   0 0 0 1;
     0 0   0    0  0 0    0     0  -1 0 0 0;
     0 0   0    0  1 0    0     0   0 0 0 0;
     L 1   0    0  0 0    0     0   0 0 0 0;
     0 0  L^3  L^2 L 1    0     0   0 0 0 0;
     0 0   0    0  0 0   L^3   L^2  L 1 0 0;
     0 0   0    0  0 0    0     0   0 0 L 1;
     0 0   0    0  0 0 -3*L^2 -2*L -1 0 0 0;
     0 0 3*L^2 2*L 1 0    0     0   0 0 0 0];

  m=inv(C)*u;

  for i=1:n
    eci(i,1)=((i-1)*L/(n-1))';
    x=eci(i,1);
    es(i,:)=([EA 0    0      0    0 0     0       0   0 0  0  0;
              0  0 -6*EIz    0    0 0     0       0   0 0  0  0;
              0  0    0      0    0 0  -6*EIy     0   0 0  0  0;
              0  0    0      0    0 0     0       0   0 0 GKv 0;
              0  0    0      0    0 0 -6*EIy*x -2*EIy 0 0  0  0;
              0  0 6*EIz*x 2*EIz  0 0     0       0   0 0  0  0;]*m + [-qx*x;
                                                                       -qy*x;
                                                                       -qz*x;
                                                                       -qw*x;
                                                                       -qz*x^2/2;
                                                                        qy*x^2/2])';
%
    edi(i,:)=([x 1  0   0  0 0  0   0  0 0 0 0;
               0 0 x^3 x^2 x 1  0   0  0 0 0 0;
               0 0  0   0  0 0 x^3 x^2 x 1 0 0;
               0 0  0   0  0 0  0   0  0 0 x 1]*m + [-qx*x^2/2/EA;
                                                      qy*x^4/24/EIz;
                                                      qz*x^4/24/EIy;
                                                     -qw*x^2/2/GKv])';
  end;
%----------------------------------end----------------------------------
