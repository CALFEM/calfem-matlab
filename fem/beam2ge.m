 function [Ke,fe]=beam2ge(ex,ey,ep,QX,eq)
% Ke=beam2ge(ex,ey,ep,QX)
% [Ke,fe]=beam2ge(ex,ey,ep,QX,eq)
%-------------------------------------------------------------
%    PURPOSE
%       Compute the element stiffness matrix for a two dimensional
%       beam element with respect to geometric nonlinearity.
%
%    INPUT: ex = [x1 x2]   element node coordinates
%           ey = [y1 y2]
%  
%           ep = [E A I]   element properties;
%                            E: Young's modulus
%                            A: Cross section area
%                            I: Moment of inertia
%
%           QX:             axial force in the beam.
% 
%           eq = [qY]       distributed transverse load
% 
%    OUTPUT: Ke:            element stiffness matrix (6 x 6)
%            fe:            element load vector (6 x 1)  
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom     2015-12-17
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

  K0le=[DEA/L   0            0      -DEA/L      0          0 ;
         0   12*DEI/L^3   6*DEI/L^2  0   -12*DEI/L^3  6*DEI/L^2;
         0   6*DEI/L^2    4*DEI/L    0   -6*DEI/L^2   2*DEI/L;
       -DEA/L  0            0       DEA/L      0          0 ;
         0   -12*DEI/L^3 -6*DEI/L^2  0   12*DEI/L^3  -6*DEI/L^2;
         0   6*DEI/L^2    2*DEI/L    0   -6*DEI/L^2   4*DEI/L];
    
  Ksle=QX/(30*L)*[0     0     0     0     0     0;
                  0    36   3*L     0   -36   3*L;
                  0   3*L 4*L^2     0  -3*L  -L^2;
                  0     0     0     0     0     0 ;
                  0   -36  -3*L     0    36  -3*L;
                  0   3*L  -L^2     0  -3*L 4*L^2];
          
  fle=qY*L*[0 1/2 L/12 0 1/2 -L/12]';

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

  Kle=K0le+Ksle; Ke=G'*Kle*G;
  fe=G'*fle;
%--------------------------end--------------------------------

 
