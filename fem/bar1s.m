 function [es,edi,eci]=bar1s(ex,ep,ed,eq,n)
% es=bar1s(ex,ep,ed)
% es=bar1s(ex,ep,ed,eq)
% [es,edi]=bar1s(ex,ep,ed,eq,n)
% [es,edi,eci]=bar1s(ex,ep,ed,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in a one dimensional bar element. 
% 
%    INPUT:  ex = [x1 x2]      element node coordinates
%
%            ep = [E A]            element properties,
%                                     E: Young's modulus
%                                     A: cross section area
%
%            ed = [u1 u2]             element displacement vector
%
%            eq = [qx]           distributed load, local direction
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

  E=ep(1); A=ep(2);
  DEA=E*A;

  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
  end
   
  qx=0;  if nargin>3;  qx=eq(1); end 
  ne=2;  if nargin>4;  ne=n; end;
     
  L=ex(2)-ex(1);
  
  C1=[1 0; -1/L 1/L];
  C1a=C1*ed';
  
  x=[0:L/(ne-1):L]';   
  zero=zeros(size(x));    one=ones(size(x));
  
  u=[one x]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u-(x.^2-L*x)*qX/(2*DEA);
    du=du-(2*x-L)*qX/(2*DEA);
  end;
  
  N=DEA*du;  
  es=N;
  edi=u;
  eci=x;
%--------------------------------- end -------------------------------------
 
