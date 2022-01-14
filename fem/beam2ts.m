 function [es,edi,eci]=beam2ts(ex,ey,ep,ed,eq,n)
% es=beam2ts(ex,ey,ep,ed)
% es=beam2ts(ex,ey,ep,ed,eq)
% [es,edi]=beam2ts(ex,ey,ep,ed,eq,n)
% [es,edi,eci]=beam2ts(ex,ey,ep,ed,eq,n)
%---------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in two dimensional Timoshenko beam 
%      element (beam2te). 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]      element node coordinates
%
%            ep = [E G A I ks] element properties,
%                              E:  Young's modulus
%                              G:  shear modulus
%                              A:  cross section area
%                              I:  moment of inertia
%                              ks: shear correction factor
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
%            edi = [ u1 v1 teta1;  element displacements, local directions,
%                    u2 v2 teta2;  in n points along the beam, dim(edi)= n x 3
%                   ............]    
%                   (Note! For Timoshenko beam element the rotation of the cross 
%                    section is not equal to dv/dx) 
%
%            eci = [x1;         evaluation points on the local x-axis 
%                   x2;      
%                   .......] 
%-------------------------------------------------------------------------
 
% LAST MODIFIED: O Dahlblom    2021-11-05
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  Gm=ep(2);   A=ep(3);  I=ep(4);  ks=ep(5);
  DEA=E*A; 
  DEI=E*I; 
  DGAks=Gm*A*ks;
  alpha=DEI/DGAks;
%  
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
  C2=1/(L^2+12*alpha)*[L^2+12*alpha 0              0          0;
                       -12*alpha/L  L^2+6*alpha    12*alpha/L -6*alpha;
                       -3           -2*L-6*alpha/L 3          -L+6*alpha/L;
                       2/L          1              -2/L       1];
  C2a=C2*a2;
     
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));
  
  u=[one X]*C1a;
  du=[zero one]*C1a;
  if DEA~=0; 
    u=u-(X.^2-L*X)*qX/(2*DEA);
    du=du-(2*X-L)*qX/(2*DEA);
  end; 
  
  v=[one X X.^2 X.^3]*C2a;
  dv=[zero one 2*X 3*X.^2]*C2a;
  theta=[zero one 2*X 3*X.^2+6*alpha]*C2a;
  dtheta=[zero zero 2*one 6*X]*C2a;
  if DEI~=0;
    v=v+(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEI)+(-X.^2/2+L*X/2)*qY/DGAks;
    dv=dv+(2*X.^3-3*L*X.^2+L^2*X)*qY/(12*DEI)+(-X+L/2)*qY/DGAks;
    theta=theta+(2*X.^3-3*L*X.^2+L^2*X)*qY/(12*DEI);
    dtheta=dtheta+(6*X.^2-6*L*X+L^2)*qY/(12*DEI);
  end;
   
  N=DEA*du; M=DEI*dtheta; V=DGAks*(dv-theta); 
  es=[N V M];
  edi=[u v theta];
  eci=X;
%--------------------------end--------------------------------
 
