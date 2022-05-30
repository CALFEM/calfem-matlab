 function [es,edi,eci]=beam2ws(ex,ey,ep,ed,eq,n)
% es=beam2ws(ex,ey,ep,ed)
% es=beam2ws(ex,ey,ep,ed,eq)
% [es,edi]=beam2ws(ex,ey,ep,ed,eq,n)
% [es,edi,eci]=beam2ws(ex,ey,ep,ed,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in a two dimensional beam element
%      on elastic foundation. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]            element node coordinates
%
%            ep = [E A I kX kY]      element properties,
%                                     E: Young's modulus
%                                     A: cross section area
%                                     I: moment of inertia
%                                     kX: axial foundation stiffness
%                                     kY: transversal found. stiffness
%
%            ed = [u1 ... u6]       element displacement vector
%
%            eq = [qX qY]           distributed loads, local directions
%          
%    OUTPUT: es = [ N1 V1 M1 ;  section forces, local directions, in 
%                   N2 V2 M2 ;  n points along the beam, dim(es)= n x 3
%                   .........]  
%           
%            edi = [ u1 v1 ;    element displacements, local directions,
%                    u2 v2 ;    in n points along the beam, dim(edi)= n x 2
%                   .......]    
%
%            eci = [x1;         evaluation points on the local x-axis 
%                   x2;      
%                   .......] 
%
% -------------------------------------------------------------------------
 
% LAST MODIFIED: O Dahlblom 2022-05-30
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  E=ep(1);  A=ep(2);  I=ep(3);  kX=ep(4);  kY=ep(5);
  DEA=E*A; DEI=E*I;
 
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qX=0; qY=0;  if nargin>4;  qX=eq(1); qY=eq(2); end 
    
  ne=2;        if nargin>5;  ne=n; end;
  
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);

  nxX=dx/L;
  nyX=dy/L;
  nxY=-dy/L;
  nyY=dx/L;
  G=[nxX  nyX   0    0    0   0;
     nxY  nyY   0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   nxX  nyX  0;
      0    0    0   nxY  nyY  0;
      0    0    0    0    0   1];

  edl=G*ed'
  
  a1=[edl(1); edl(4)];
  C1=[1 0;-1/L 1/L];
  C1a=C1*a1;
  
  a2=[edl(2); edl(3); edl(5); edl(6)];
  C2=[1       0    0       0;
      0       1    0       0;
     -3/(L^2) -2/L 3/(L^2) -1/L;
     2/(L^3) 1/(L^2) -2/(L^3) 1/(L^2)];
  C2a=C2*a2;
   
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));
  
  u=[one X]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u+kX/DEA*[(X.^2-L*X)/2 (X.^3-L^2*X)/6]*C1a-(X.^2-L*X)*qX/(2*DEA);
    du=du+kX/DEA*[(2*X-L)/2 (3*X.^2-L^2)/6]*C1a-(2*X-L)*qX/(2*DEA);
  end; 
  v=[one X X.^2 X.^3]*C2a;
  d2v=[zero zero 2*one 6*X]*C2a;
  d3v=[zero zero zero 6*one]*C2a;
  if DEI~=0;
    v=v-kY/DEI*[(X.^4-2*L*X.^3+L^2*X.^2)/24 (X.^5-3*L^2*X.^3+2*L^3*X.^2)/120 (X.^6-4*L^3*X.^3+3*L^4*X.^2)/360 (X.^7-5*L^4*X.^3+4*L^5*X.^2)/840]*C2a...
    +(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEI);
    d2v=d2v-kY/DEI*[(6*X.^2-6*L*X+L^2)/12 (10*X.^3-9*L^2*X+2*L^3)/60 (5*X.^4-4*L^3*X+L^4)/60 (21*X.^5-15*L^4*X+4*L^5)/420]*C2a...
    +(6*X.^2-6*L*X+L^2)*qY/(12*DEI);
    d3v=d3v-kY/DEI*[(2*X-L)/2 (10*X.^2-3*L^2)/20 (5*X.^3-L^3)/15 (7*X.^4-L^4)/28]*C2a...
    +(2*X-L)*qY/(2*DEI);
  end;   
  N=DEA*du; M=DEI*d2v; V=-DEI*d3v; 
  es=[N V M];
  edi=[u v];
  eci=X;
%--------------------------------- end -------------------------------------
 
