function [es,edi,eci]=bar3s(ex,ey,ez,ep,ed,eq,n)
% es=bar3s(ex,ey,ez,ep,ed)
% es=bar3s(ex,ey,ez,ep,ed,eq)
% [es,edi]=bar3s(ex,ey,ez,ep,ed,eq,n)
% [es,edi,eci]=bar2s(ex,ey,ep,ed,eq,n)
%-------------------------------------------------------------
% PURPOSE
%  Compute normal force in three dimensional bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]       
%         ez = [z1 z2]          element node coordinates
%    
%         ep = [E A]            element properties   
%                                  E : Young's modulus
%                                  A : Cross section area
%
%         ed : [u1 ... u6]      element displacements
%
%         eq = [qX]             distributed load, local direction
%
%
%    OUTPUT: es = [N1;
%                  N2 ]             section forces, local directions
%
%            edi = [ u1 ;           element displacements, local directions,
%                    u2 ;           in n points along the bar, dim(es)= n x 1
%                   ...]    
%
%            eci = [ x1  ;          local x-coordinates of the evaluation 
%                    x2 ;           points, (x1=0 and xn=L)
%                    ...]
%
%-------------------------------------------------------------

% LAST MODIFIED: O. Dahlblom    2021-09-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  
  E=ep(1);  A=ep(2);  
  DEA=E*A;
 
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qX=0;  if nargin>5;  qX=eq(1); end 
    
  ne=2;  if nargin>6;  ne=n; end;
  
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  dy=ez(2)-ez(1);
  L=sqrt(dx*dx+dy*dy+dz*dz);
% b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
% L=sqrt(b'*b);

  nxX=dx/L;
  nyX=dy/L;
  nzX=dz/L;
  G=[nxX nyX nzX  0   0  0 ;  
      0   0   0 nxX nyX nzX];

  a1=G*ed';
  
  C1=[1 0;-1/L 1/L];
  C1a=C1*a1;
  
  X=[0:L/(ne-1):L]';
  zero=zeros(size(X));    one=ones(size(X));
  
  u=[one X]*C1a-(X.^2-L*X)*qX/(2*DEA);
  du=[zero one]*C1a-(2*X-L)*qX/(2*DEA);
   
  N=DEA*du;
  es=N;
  edi=u;
  eci=X;
%--------------------------end--------------------------------
