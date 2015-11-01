 function [es,edi,eci]=beam2s(ex,ey,ep,ed,eq,n)
% es=beam2s(ex,ey,ep,ed)
% es=beam2s(ex,ey,ep,ed,eq)
% [es,edi,eci]=beam2s(ex,ey,ep,ed,eq,n)
%---------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in two dimensional beam element (beam2e). 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]     element node coordinates
%
%            ep = [E A I]     element properties,
%                              E:  Young's modulus
%                              A:  cross section area
%                              I:  moment of inertia
%
%            ed = [u1 ... u6] element displacements
%
%            eq = [qx qy]     distributed loads, local directions 
%
%            n : number of evaluation points ( default=2 )
%          
%    OUTPUT: es = [ N1 V1 M1 ;  section forces, local directions, in 
%                   N2 V2 M2 ;  n points along the beam, dim(es)= n x 3
%                   .........]  
%           
%            edi = [ u1 v1 ;    element displacements, local directions,
%                    u2 v2 ;    in n points along the beam, dim(es)= n x 2
%                   .......]    
%
%            eci = [ x1  ;      local x-coordinates of the evaluation 
%                    x2 ;       points, (x1=0 and xn=L)
%                    ...]
%-------------------------------------------------------------------------
 
% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  EA=ep(1)*ep(2); EI=ep(1)*ep(3);
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);

  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qx=0; qy=0;  if nargin>4;  qx=eq(1); qy=eq(2); end 
    
  ne=2;        if nargin>5;  ne=n; end;
     
  C=[0   0   0    1   0   0;
     0   0   0    0   0   1;
     0   0   0    0   1   0;
     L   0   0    1   0   0;
     0   L^3  L^2 0   L   1;
     0 3*L^2 2*L  0   1   0];

  n=b/L;

  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
 
  M=inv(C)*(G*ed'-[0 0 0 -qx*L^2/(2*EA) qy*L^4/(24*EI) qy*L^3/(6*EI)]' );
  
  A=[M(1) M(4)]';  B=[M(2) M(3) M(5) M(6)]';
  
  x=[0:L/(ne-1):L]';   zero=zeros(size(x));    one=ones(size(x));
  
  u=[x one]*A-(x.^2)*qx/(2*EA);
  du=[one zero]*A-x*qx/EA;
  v=[x.^3 x.^2 x one]*B+(x.^4)*qy/(24*EI); 
% dv=[3*x.^2 2*x one zero]*B+(x.^3)*qy/(6*EI);
  d2v=[6*x 2*one zero zero]*B+(x.^2)*qy/(2*EI);
  d3v=[6*one zero zero zero]*B+x*qy/EI;
  
  N=EA*du; M=EI*d2v; V=-EI*d3v; 
  es=[N V M];
  edi=[u v];
  eci=x;
%--------------------------end--------------------------------
 
