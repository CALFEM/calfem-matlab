 function [Ke,fe]=beam2gxe(ex,ey,ep,QX,eq)
% Ke=beam2gxe(ex,ey,ep,QX)
% [Ke,fe]=beam2gxe(ex,ey,ep,QX,eq)
%-------------------------------------------------------------
%    PURPOSE
%       Compute the element stiffness matrix for a two dimensional
%       beam element with respect to geometric nonlinearity with exact solution.
%
%    INPUT: ex = [x1 x2]   element node coordinates
%           ey = [y1 y2]
%  
%           ep = [E A I]   element properties;
%                            E: Young's modulus
%                            A: Cross section area
%                            I: Moment of inertia
%
%           QX:  		axial force in the beam.
% 
%           eq = [qY]		distributed transverse load
%              
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%            fe : element load vector (6 x 1)  
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2021-06-21
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  A=ep(2);  I=ep(3); 
  DEA=E*A; DEI=E*I;

  qY=0;  
  if nargin>4; 
     if length(eq)>1
        disp('eq should be a scalar!!!')
        return 
     end
     qY=eq(1); 
  end

  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);
 
  eps=1e-12;
  if QX<-eps*DEI/L^2
     kL=sqrt(-QX/DEI)*L;
         
     f1=(kL/2)/tan(kL/2);
     f2=kL^2/(12*(1-f1));
     f3=f1/4+3*f2/4;
     f4=-f1/2+3*f2/2;
     f5=f1*f2;
     
     h=6*(2/kL^2-(1+cos(kL))/(kL*sin(kL)));
  elseif QX>eps*DEI/L^2
     kL=sqrt(QX/DEI)*L;
         
     f1=(kL/2)/tanh(kL/2);
     f2=-kL^2/(12*(1-f1));
     f3=f1/4+3*f2/4;
     f4=-f1/2+3*f2/2;
     f5=f1*f2;

     h=-6*(2/kL^2-(1+cosh(kL))/(kL*sinh(kL)));
  else
     f1=1; f2=1; f3=1; f4=1; f5=1; h=1;
  end

  Kle=[DEA/L  0            0        -DEA/L      0          0 ;
       0  12*DEI*f5/L^3   6*DEI*f2/L^2  0 -12*DEI*f5/L^3 6*DEI*f2/L^2;
       0  6*DEI*f2/L^2    4*DEI*f3/L    0  -6*DEI*f2/L^2   2*DEI*f4/L;
       -DEA/L  0            0         DEA/L      0          0 ;
       0  -12*DEI*f5/L^3 -6*DEI*f2/L^2  0  12*DEI*f5/L^3 -6*DEI*f2/L^2;
       0  6*DEI*f2/L^2    2*DEI*f4/L    0  -6*DEI*f2/L^2   4*DEI*f3/L];

  fle=qY*L*[0 1/2 L*h/12 0 1/2 -L*h/12]';

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

  Ke=G'*Kle*G;  
  fe=G'*fle;
%--------------------------end--------------------------------

 
