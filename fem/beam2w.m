 function [Ke,fe]=beam2w(ex,ey,ep,eq)
% Ke=beam2w(ex,ey,ep)
% [Ke,fe]=beam2w(ex,ey,ep,eq)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute the stiffness matrix for a two dimensional beam element 
%      on elastic foundation.
%
%    INPUT:    ex = [x1 x2]
%              ey = [y1 y2]         element node coordinates
% 
%              ep = [E A I ka kt]   element properties;
%                                     E: Young's modulus
%                                     A: cross section area
%                                     I: moment of inertia
%                                     ka: axial foundation stiffness
%                                     kt: transversal found. stiffness
%
%              eq = [qx qy]           distributed loads (local directions)
%            
%    OUTPUT:   Ke : beam stiffness matrix (6 x 6)
%
%              fe : element load vector (6 x 1)
%-------------------------------------------------------------------------
  
% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%--------------------------------------------------------------------------

  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;
%
  E=ep(1);  A=ep(2);  I=ep(3);  ka=ep(4);  kt=ep(5);
% 
  qx=0; qy=0;  if nargin>3; qx=eq(1); qy=eq(2); end
%
  K1 =[E*A/L   0            0      -E*A/L      0          0 ;
         0   12*E*I/L^3   6*E*I/L^2  0   -12*E*I/L^3  6*E*I/L^2;
         0   6*E*I/L^2    4*E*I/L    0   -6*E*I/L^2   2*E*I/L;
       -E*A/L  0            0       E*A/L      0          0 ;
         0   -12*E*I/L^3 -6*E*I/L^2  0   12*E*I/L^3  -6*E*I/L^2;
         0   6*E*I/L^2    2*E*I/L    0   -6*E*I/L^2   4*E*I/L];
%         
  K2=L/420*[140*ka   0       0      70*ka    0       0     ;
             0    156*kt   22*kt*L    0    54*kt  -13*kt*L ;
             0    22*kt*L  4*kt*L^2   0   13*kt*L -3*kt*L^2;
            70*ka    0       0     140*ka    0       0     ;
             0    54*kt    13*kt*L    0   156*kt  -22*kt*L ;
             0   -13*kt*L -3*kt*L^2   0  -22*kt*L  4*kt*L^2];
%   
  Kle=K1+K2;
%   
  fle=L*[qx/2 qy/2 qy*L/12 qx/2 qy/2 -qy*L/12]';
%
  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
%
  Ke=G'*Kle*G;   fe=G'*fle;    
%------------------------------------end-------------------------------------
  
