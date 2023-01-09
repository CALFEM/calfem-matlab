 function [es,edi,eci]=beam1ws(ex,ep,ed,eq,n)
% es=beam1ws(ex,ep,ed)
% es=beam1ws(ex,ep,ed,eq)
% [es,edi]=beam1ws(ex,ep,ed,eq,n)
% [es,edi,eci]=beam1ws(ex,ep,ed,eq,n)
%---------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in one dimensional beam element 
%      on elastic foundation (beam1we). 
% 
%    INPUT:  ex = [x1 x2]     element node coordinates
%
%            ep = [E I kY]     element properties,
%                              E:  Young's modulus
%                              I:  moment of inertia
%                              kY: transversal found. stiffness
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
 
% LAST MODIFIED: O Dahlblom 2021-09-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  I=ep(2);  kY=ep(3);
  DEI=E*I;
  
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qY=0;  if nargin>3;  qY=eq(1); end 
    
  ne=2;        if nargin>4;  ne=n; end;
     
  dx=ex(2)-ex(1);
  L=abs(dx);
   
  a1=ed';
  C2=[1       0    0       0;
      0       1    0       0;
     -3/(L^2) -2/L 3/(L^2) -1/L;
     2/(L^3) 1/(L^2) -2/(L^3) 1/(L^2)];
  C2a=C2*a1;
     
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));
  
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

  M=DEI*d2v; V=-DEI*d3v; 
  es=[V M];
  edi=[v];
  eci=X;
%--------------------------end--------------------------------
 
