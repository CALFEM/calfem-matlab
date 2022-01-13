 function [es,edi,eci]=bar1ws(ex,ep,ed,eq,n)
% es=bar1ws(ex,ep,ed)
% es=bar1ws(ex,ep,ed,eq)
% [es,edi]=bar1ws(ex,ep,ed,eq,n)
% [es,edi,eci]=bar1ws(ex,ep,ed,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in a one dimensional bar element
%      with elastic support. 
% 
%    INPUT:  ex = [x1 x2]      element node coordinates
%
%            ep = [E A ka]     element properties,
%                                E: Young's modulus
%                                A: cross section area
%                                ka: axial foundation stiffness
%
%            ed = [u1 u2]      element displacement vector
%
%            eq = [qx]         distributed load, local direction
%
%            n : number of evaluation points ( default=2 )
%          
%    OUTPUT: es = [N1;
%                  N2]       element forces, local directions 
%
%            edi = [ u1 ;    element displacements, local directions,
%                    u2 ;    in n points along the bar, dim(es)= n x 1
%                   .......]    
%
%            eci = [ x1  ;      local x-coordinates of the evaluation 
%                    x2 ;       points, (x1=0 and xn=L)
%                    ...]
% -------------------------------------------------------------------------
 
% LAST MODIFIED: O Dahlblom  2021-02-25
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
  end
  
  L=ex(2)-ex(1);
  
  qx=0;  if nargin>3;  qx=eq(1); end 
    
  ne=2;  if nargin>4;  ne=n; end;  

  E=ep(1); A=ep(2); ka=ep(3);
  DEA=E*A;
  
  C1=[1 0;-1/L 1/L];
  C1a=C1*ed';
  
  x=[0:L/(ne-1):L]';   zero=zeros(size(x));    one=ones(size(x));
  
  u=[one x]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u+ka/DEA*[(x.^2-L*x)/2 (x.^3-L^2*x)/6]*C1a-(x.^2-L*x)*qx/(2*DEA);
    du=du+ka/DEA*[(2*x-L)/2 (3*x.^2-L^2)/6]*C1a-(2*x-L)*qx/(2*DEA);
  end; 

  N=DEA*du;  
  es=[N];
  edi=[u];
  eci=x;
%--------------------------------- end -------------------------------------
 
