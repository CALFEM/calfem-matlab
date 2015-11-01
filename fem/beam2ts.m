function [es,edi,eci]=beam2ts(ex,ey,ep,ed,eq,n)
%        [es,edi,eci]=beam2ts(ex,ey,ep,ed,eq,n) 
%-------------------------------------------------------------
%    PURPOSE
%      Compute section forces in two dimensional
%      Timoshenko beam element (beam2te). 
% 
%    INPUT:  ex = [x1 x2];        
%            ey = [y1 y2];        element node coordinates
%
%            ep = [E G A I ks ];  element properties,
%                                   E: Young's modulus
%                                   G: shear modulus
%                                   A: cross section area
%                                   I: moment of inertia
%                                  ks: shear correction factor
%
%            ed = [u1 ... u6] element displacements
%                 
%
%            eq = [qx qy]     distributed loads, local directions
%
%            n : number of evaluation points (default=2)
%
%    OUTPUT: es = [ N1 V1 M1 ;  section forces, local directions, in
%                   N2 V2 M2 ;  n points along the beam, dim(es)= n x 3
%                   .........]
%
%            edi = [ u1 v1 teta1;    element displacements, local directions,
%                    u2 v2 teta2;    and rotation of cross section at
%                   ............]    n points along the beam, dim(es)= n x 3
%
%(Note! Rotation of the cross section is not equal to dv/dx for Timoshenko beam element)
%                                
%            eci = [ x1  ;      local x-coordinates of the evaluation
%                    x2 ;       points, (x1=0 and xn=L)
%                    ...]
%-------------------------------------------------------------------------

% LAST MODIFIED: E Serrano   1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  if nargin==5; n=2; end;
  if nargin==4; n=2; eq=[0 0]; end;
  ne=n;
%
  EA=ep(1)*ep(3);
  EI=ep(1)*ep(4);
  GAK=ep(2)*ep(3)*ep(5);
  alfa=EI/GAK;
%
  b=[ex(2)-ex(1);ey(2)-ey(1)];
  L=sqrt(b'*b); n=b/L;
%
  qx=eq(1);qy=eq(2);
%
  C=[0       0          0    1   0   0;
     0       0          0    0   0   1;
     0     6*alfa       0    0   1   0;
     L       0          0    1   0   0;
     0       L^3       L^2   0   L   1;
     0  3*(L^2+2*alfa) 2*L   0   1   0];
%
  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
%
  M=inv(C)*(G*ed'-[0 0 0 -qx*L^2/(2*EA) qy*L^4/(24*EI)-qy*L^2/(2*GAK) qy*L^3/(6*EI)]' );
  C2=[M(1) M(4)]';
  C4=[M(2) M(3) M(5) M(6)]';
%
  x(1:ne,1)=[0:(L/(ne-1)):L]';
%
  one=ones(size(x));
  zero=zeros(size(x));
%
  u=[x one]*C2-qx/(2*EA)*x.^2;
  du=[one zero]*C2-qx*x/EA;
%
  v=[x.^3 x.^2 x one]*C4+qy/(24*EI)*(x.^4)-qy/(2*GAK)*(x.^2);
  dv=[3*x.^2 2*x one zero]*C4+qy*(x.^3)/(6*EI)-qy*x/GAK;
%
  teta=[3*(x.^2+2*alfa*one) 2*x one zero]*C4+qy*(x.^3)/(6*EI);
  dteta=[6*x 2*one zero zero]*C4+qy*(x.^2)/(2*EI); 
%
  N=EA*du;
  M=EI*dteta;
  V=GAK*(dv-teta);
%
  es=[N V M];
  edi=[u v teta];
  eci=[x];
%--------------------------- end -----------------------------
