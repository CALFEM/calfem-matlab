 function [es,edi,eci]=beam1s(ex,ep,ed,eq,n)
% es=beam1s(ex,ep,ed)
% es=beam1s(ex,ep,ed,eq)
% [es,edi]=beam1s(ex,ep,ed,eq,n)
% [es,edi,eci]=beam1s(ex,ep,ed,eq,n)
%---------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in two dimensional beam element (beam2e). 
% 
%    INPUT:  ex = [x1 x2]     element node coordinates
%
%            ep = [E I]     element properties,
%                              E:  Young's modulus
%                              I:  moment of inertia
%
%            ed = [u1 ... u4] element displacements
%
%            eq = [qy]     distributed loads, local directions 
%
%            n : number of evaluation points ( default=2 )
%          
%    OUTPUT: es = [V1 M1 ;  section forces, local directions, in 
%                  V2 M2 ;  n points along the beam, dim(es)= n x 2
%                   .........]  
%           
%            edi = [v1 ;    element displacements, local directions,
%                   v2 ;    in n points along the beam, dim(edi)= n x 1
%                   .......]
%
%            eci = [x1;     evaluation points on the local x-axis 
%                   x2;      
%                   .......] 
%-------------------------------------------------------------------------
 
% LAST MODIFIED: O Dahlblom    2021-09-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  I=ep(2); 
  DEI=E*I;
  
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qY=0;  if nargin>3;  qY=eq(1); end 
    
  ne=2;        if nargin>4;  ne=n; end;
     
  L=ex(2)-ex(1);
     
  a1=ed';
  C2=[1       0    0       0;
      0       1    0       0;
     -3/(L^2) -2/L 3/(L^2) -1/L;
     2/(L^3) 1/(L^2) -2/(L^3) 1/(L^2)];
  C2a=C2*a1;
     
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));
  
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

  M=DEI*d2v; V=-DEI*d3v; 
  es=[V M];
  edi=[v];
  eci=X;
%--------------------------end--------------------------------
 
