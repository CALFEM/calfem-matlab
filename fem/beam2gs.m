 function [es,QX,edi,eci]=beam2gs(ex,ey,ep,ed,QX,eq,n)
% [es,QX]=beam2gs(ex,ey,ep,ed,QX)
% [es,QX]=beam2gs(ex,ey,ep,ed,QX,eq)
% [es,QX,edi]=beam2gs(ex,ey,ep,ed,QX,eq,n)
% [es,QX,edi,eci]=beam2gs(ex,ey,ep,ed,QX,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Calculate section forces in a two dimensional nonlinear
%      beam element (beam2ge).
% 
%    INPUT:  ex = [x1 x2]           element node coordinates 
%            ey = [y1 y2]           
%
%            ep = [E A I]           element properties;
%                                     E: Young's modulus
%                                     A: cross section area
%                                     I: moment of inertia
%
%            ed = [u1 ... u6]       element displacement vector
%
%            QX:                    axial force in the beam
%
%            eq = [qY]              distributed transverse load
%
%            n:                     number of evaluation points ( default=2 )
%          
%    OUTPUT: es = [N1 V1 M1 ;
%                  N2 V2 M2 ]       element forces, local directions
%
%            QX:                    axial force
%
%            edi = [ u1 v1 ;        element displacements, local directions,
%                    u2 v2 ;        in n points along the beam, dim(es)= n x 2
%                   .......]    
%
%            eci = [x1;             evaluation points on the local x-axis 
%                   x2;      
%                   .......] 
%
% -------------------------------------------------------------------------
 
% LAST MODIFIED: O. Dahlblom    2021-09-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  E=ep(1);  A=ep(2);  I=ep(3);  
  DEA=E*A; DEI=E*I;
 
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
   end
  
  qY=0;  
  if nargin>5; 
     if length(eq)>1
        disp('eq should be a scalar!!!')
        return 
     end
     qY=eq(1); 
  end
    
  ne=2;  if nargin>6;  ne=n; end;
  
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
  v=[one X X.^2 X.^3]*C2a;  
  dv=[zero one 2*X 3*X.^2]*C2a;
  d2v=[zero zero 2*one 6*X]*C2a;
  d3v=[zero zero zero 6*one]*C2a;
  if DEI~=0;
    v=v+QX/DEI*[zero zero (X.^4-2*L*X.^3+L^2*X.^2)/12 (X.^5-3*L^2*X.^3+2*L^3*X.^2)/20]*C2a...
    +(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEI);
    dv=dv+QX/DEI*[zero zero (2*X.^3-3*L*X.^2+L^2*X)/6 (5*X.^4-9*L^2*X.^2+4*L^3*X)/20]*C2a...
    +(2*X.^3-3*L*X.^2+L^2*X)*qY/(12*DEI);
    d2v=d2v+QX/DEI*[zero zero (6*X.^2-6*L*X+L^2*one)/6 (10*X.^3-9*L^2*X+2*L^3*one)/10]*C2a...
    +(6*X.^2-6*L*X+L^2)*qY/(12*DEI);
    d3v=d3v+QX/DEI*[zero zero (2*X-L*one) (30*X.^2-9*L^2*one)/10]*C2a...
    +(2*X-L)*qY/(2*DEI);
  end;
   
  QX=DEA*du(1);
  M=DEI*d2v; V=-DEI*d3v;
  N=QX+dv.*V;
  es=[N V M];
  edi=[u v];
  eci=X;
%--------------------------------- end -------------------------------------
 
