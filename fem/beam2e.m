function [Ke,fe]=beam2e(ex,ey,ep,eq);
% Ke=beam2e(ex,ey,ep)
% [Ke,fe]=beam2e(ex,ey,ep,eq)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the stiffness matrix for a two dimensional beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%            ep = [E A I]       element properties
%                                  E: Young's modulus
%                                  A: Cross section area
%                                  I: Moment of inertia
%            eq = [qX qY]       distributed loads, local directions
% 
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%            fe : element load vector (6 x 1)
%--------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2015-08-07
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  A=ep(2);  I=ep(3); 
  DEA=E*A; DEI=E*I;
  
  qX=0; qY=0;  if nargin>3; qX=eq(1); qY=eq(2); end

  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);
  
  Kle=[DEA/L   0            0      -DEA/L      0          0 ;
         0   12*DEI/L^3   6*DEI/L^2  0   -12*DEI/L^3  6*DEI/L^2;
         0   6*DEI/L^2    4*DEI/L    0   -6*DEI/L^2   2*DEI/L;
       -DEA/L  0            0       DEA/L      0          0 ;
         0   -12*DEI/L^3 -6*DEI/L^2  0   12*DEI/L^3  -6*DEI/L^2;
         0   6*DEI/L^2    2*DEI/L    0   -6*DEI/L^2   4*DEI/L];
   
  fle=L*[qX/2 qY/2 qY*L/12 qX/2 qY/2 -qY*L/12]';

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

  Ke=G'*Kle*G;   fe=G'*fle; 
%--------------------------end--------------------------------
