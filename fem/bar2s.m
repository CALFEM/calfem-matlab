 function [es,edi,eci]=bar2s(ex,ey,ep,ed,eq,n)
% es=bar2s(ex,ey,ep,ed)
% es=bar2s(ex,ey,ep,ed,eq)
% [es,edi]=bar2s(ex,ey,ep,ed,eq,n)
% [es,edi,eci]=bar2s(ex,ey,ep,ed,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Calculate section forces in a two dimensional bar element (bar2e).
% 
%    INPUT:  ex = [x1 x2]           element node coordinates 
%            ey = [y1 y2]           
%
%            ep = [E A]           element properties;
%                                     E: Young's modulus
%                                     A: cross section area
%
%            ed = [u1 ... u4]       element displacement vector
%
%            eq = [qX]              distributed load, local direction
%
%            n:                     number of evaluation points ( default=2 )
%          
%    OUTPUT: es = [N1 ;             section forces, local directions
%                  N2 ;
%                  ...]              
%
%            edi = [ u1 ;           element displacements, local directions,
%                    u2 ;           in n points along the bar, dim(es)= n x 1
%                   ...]    
%
%            eci = [ x1 ;      local x-coordinates of the evaluation 
%                    x2 ;       points, (x1=0 and xn=L)
%                    ...]
% -------------------------------------------------------------------------
 
% LAST MODIFIED: O. Dahlblom    2015-12-04
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  
  E=ep(1);  A=ep(2);  
  DEA=E*A;
 
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qX=0;  if nargin>4;  qX=eq(1); end; 
    
  ne=2;  if nargin>5;  ne=n; end;
  
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);
   
  nxX=dx/L;
  nyX=dy/L;
  G=[nxX  nyX   0    0;
      0    0    nxX  nyX];
 
  a1=G*ed';
  
  C1=[1 0;-1/L 1/L];
  C1a=C1*a1;
  
  X=[0:L/(ne-1):L]';
  zero=zeros(size(X));    one=ones(size(X));
  
  u=[one X]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u-(X.^2-L*X)*qX/(2*DEA);
    du=du-(2*X-L)*qX/(2*DEA);
  end; 
  
  N=DEA*du;
  es=N;
  edi=u;
  eci=X;
%--------------------------------- end -------------------------------------
 
