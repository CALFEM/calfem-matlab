 function [es,edi,eci]=beam2s(ex,ey,ep,ed,eq,n)
% es=beam2s(ex,ey,ep,ed)
% es=beam2s(ex,ey,ep,ed,eq)
% [es,edi]=beam2s(ex,ey,ep,ed,eq,n)
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
%            eq = [qX qY]     distributed loads, local directions 
%
%            n : number of evaluation points ( default=2 )
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
%-------------------------------------------------------------------------
 
% LAST MODIFIED: O Dahlblom    2021-09-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  A=ep(2);  I=ep(3); 
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
 
  edl=G*ed';
  
  a1=[edl(1); edl(4)];
  C1=[ 1   0;
      -1/L 1/L];
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
    u=u-(X.^2-L*X)*qX/(2*DEA);
    du=du-(2*X-L)*qX/(2*DEA);
  end; 
  
  v=[one X X.^2 X.^3]*C2a;
% dv=[zero one 2*X 3*X.^2]*C2a;
  d2v=[zero zero 2*one 6*X]*C2a;
  d3v=[zero zero zero 6*one]*C2a;
  if DEI~=0;
    v=v+(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEI);
%   dv=dv+(2*X.^3-3*L*X.^2+L^2*X)*qY/(12*DEI);
    d2v=d2v+(6*X.^2-6*L*X+L^2)*qY/(12*DEI);
    d3v=d3v+(2*X-L)*qY/(2*DEI);
  end;
   
  N=DEA*du; M=DEI*d2v; V=-DEI*d3v; 
  es=[N V M];
  edi=[u v];
  eci=X;
%--------------------------end--------------------------------
 
