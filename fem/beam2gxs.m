 function [es,QX,edi,eci]=beam2gxs(ex,ey,ep,ed,QX,eq,n)
% [es,QX]=beam2gxs(ex,ey,ep,ed,QX)
% [es,QX]=beam2gxs(ex,ey,ep,ed,QX,eq)
% [es,QX,edi]=beam2gxs(ex,ey,ep,ed,QX,eq,n)
% [es,QX,edi,eci]=beam2gxs(ex,ey,ep,ed,QX,eq,n)
%-------------------------------------------------------------------------
%    PURPOSE
%      Calculate section forces in a two dimensional nonlinear
%      beam element (beam2gxe).
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
 
% LAST MODIFIED: O. Dahlblom    2021-09-17
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
  
  X=[0:L/(ne-1):L]';   zero=zeros(size(X));    one=ones(size(X));
  
  a1=[edl(1); edl(4)];
  C1=[ 1   0;
      -1/L 1/L];
  C1a=C1*a1;
  u=[one X]*C1a;
  du=[zero one]*C1a;
  
  a2=[edl(2); edl(3); edl(5); edl(6)];
  
  eps=1e-12;
  if QX<-eps*DEI/L^2
     k=sqrt(-QX/DEI);
     kL=k*L;
     C2=1/(k*(-2*(1-cos(kL))+kL*sin(kL)))*...
     [k*(kL*sin(kL)+cos(kL)-1) -kL*cos(kL)+sin(kL)   -k*(1-cos(kL)) -sin(kL)+kL;
     -(k^2)*sin(kL)            -k*(1-cos(kL))        (k^2)*sin(kL)  -k*(1-cos(kL));
     -k*(1-cos(kL))            kL*cos(kL)-sin(kL)     k*(1-cos(kL))  sin(kL)-kL;
     k*sin(kL)                 (kL*sin(kL)+cos(kL)-1) -k*sin(kL)    (1-cos(kL))];
     C2a=C2*a2;
     v=[one X cos(k.*X) sin(k.*X)]*C2a;  
     dv=[zero one -k*sin(k.*X) k*cos(k.*X)]*C2a;
     d2v=[zero zero -k^2*cos(k.*X) -k^2*sin(k.*X)]*C2a;
     d3v=[zero zero k^3*sin(k.*X) -k^3*cos(k.*X)]*C2a;
      if DEI~=0;
         v=v+qY*L^4/(2*DEI)*((1+cos(kL))/(kL^3*sin(kL))*(-1+cos(k.*X))+sin(k.*X)/kL^3+X.*(-1+X./L)/(kL^2*L));
         dv=dv+qY*L^3/(2*DEI)*((1+cos(kL))/(kL^2*sin(kL))*(-sin(k.*X))+cos(k.*X)/kL^2+(-1+2*X./L)/kL^2);
         d2v=d2v+qY*L^2/(2*DEI)*((1+cos(kL))/(kL*sin(kL))*(-cos(k.*X))-sin(k.*X)/kL+2/(kL^2)); 
         d3v=d3v+qY*L/(2*DEI)*((1+cos(kL))/sin(kL)*(sin(k.*X))-cos(k.*X));
      end;
  elseif QX>eps*DEI/L^2
     k=sqrt(QX/DEI);
     kL=k*L;
     C2=1/(k*(-2*(1-cosh(kL))-kL*sinh(kL)))*...
     [k*(-kL*sinh(kL)+cosh(kL)-1) -kL*cosh(kL)+sinh(kL)     -k*(1-cosh(kL)) -sinh(kL)+kL;
     (k^2)*sinh(kL)               -k*(1-cosh(kL))           -(k^2)*sinh(kL) -k*(1-cosh(kL));
     -k*(1-cosh(kL))              kL*cosh(kL)-sinh(kL)      k*(1-cosh(kL))  sinh(kL)-kL;
     -k*sinh(kL)                  (-kL*sinh(kL)+cosh(kL)-1) k*sinh(kL)      (1-cosh(kL))];
     C2a=C2*a2;
     v=[one X cosh(k.*X) sinh(k.*X)]*C2a;  
     dv=[zero one k*sinh(k.*X) k*cosh(k.*X)]*C2a;
     d2v=[zero zero k^2*cosh(k.*X) k^2*sinh(k.*X)]*C2a;
     d3v=[zero zero k^3*sinh(k.*X) k^3*cosh(k.*X)]*C2a;
      if DEI~=0;
         v=v+qY*L^4/(2*DEI)*((1+cosh(kL))/(kL^3*sinh(kL))*(-1+cosh(k.*X))-sinh(k.*X)/kL^3+X.*(1-X./L)/(kL^2*L));
         dv=dv+qY*L^3/(2*DEI)*((1+cosh(kL))/(kL^2*sinh(kL))*(sinh(k.*X))-cosh(k.*X)/kL^2+(1-2*X./L)/kL^2);
         d2v=d2v+qY*L^2/(2*DEI)*((1+cosh(kL))/(kL*sinh(kL))*(cosh(k.*X))-sinh(k.*X)/kL-2/(kL^2)); 
         d3v=d3v+qY*L/(2*DEI)*((1+cosh(kL))/sinh(kL)*(sinh(k.*X))-cosh(k.*X));
      end;
  else
     C2=[1       0       0          0;
         0       1       0          0;
        -3/(L^2) -2/L    3/(L^2)    -1/L;
        2/(L^3)  1/(L^2) -2/(L^3)   1/(L^2)];
     C2a=C2*a2;
     v=[one X X.^2 X.^3]*C2a;  
     dv=[zero one 2*X 3*X.^2]*C2a;
     d2v=[zero zero 2*one 6*X]*C2a;
     d3v=[zero zero zero 6*one]*C2a;
     if DEI~=0;
        v=v+(X.^4-2*L*X.^3+L^2*X.^2)*qY/(24*DEI);
        dv=dv+(2*X.^3-3*L*X.^2+L^2*X)*qY/(12*DEI);
        d2v=d2v+(6*X.^2-6*L*X+L^2)*qY/(12*DEI);
        d3v=d3v+(2*X-L)*qY/(2*DEI);
     end;
  end
    
  QX=DEA*du(1);
  M=DEI*d2v; V=-DEI*d3v;
  N=QX+dv.*V;
  es=[N V M];
  edi=[u v];
  eci=X;
  %--------------------------------- end -------------------------------------
 
