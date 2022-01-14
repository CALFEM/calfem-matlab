function [Ke,fe]=beam2te(ex,ey,ep,eq);
% Ke=beam2te(ex,ey,ep)
% [Ke,fe]=beam2te(ex,ey,ep,eq)
%---------------------------------------------------------------------
%    PURPOSE
%     Compute the stiffness matrix for a two dimensional Timoshenko 
%     beam element. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]       element node coordinates
%            ep = [E G A I ks   element properties
%                                  E: Young's modulus
%                                  G: shear modulus
%                                  A: Cross section area
%                                  I: Moment of inertia
%                                  ks: shear correction factor
%            eq = [qX qY]       distributed loads, local directions
% 
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%            fe : element load vector (6 x 1)
%--------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2021-11-05
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  Gm=ep(2);   A=ep(3);  I=ep(4);  ks=ep(5);
  DEA=E*A; DEI=E*I; DGAks=Gm*A*ks;
  
  qX=0; qY=0;  if nargin>3; qX=eq(1); qY=eq(2); end

  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);
  m=(12*DEI)/(L^2*DGAks);
  f1=1/(1+m);
  f2=f1*(1+m/4);
  f3=f1*(1-m/2);
  
  Kle=[DEA/L  0              0             -DEA/L  0             0 ;
         0    12*DEI*f1/L^3  6*DEI*f1/L^2  0      -12*DEI*f1/L^3 6*DEI*f1/L^2;
         0    6*DEI*f1/L^2   4*DEI*f2/L    0      -6*DEI*f1/L^2  2*DEI*f3/L;
       -DEA/L 0              0             DEA/L      0          0;
         0    -12*DEI*f1/L^3 -6*DEI*f1/L^2 0      12*DEI*f1/L^3  -6*DEI*f1/L^2;
         0    6*DEI*f1/L^2   2*DEI*f3/L    0      -6*DEI*f1/L^2  4*DEI*f2/L];
   
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
