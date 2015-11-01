function [Ke]=bar2e(ex,ey,ep)
% Ke=bar2e(ex,ey,ep)
%----------------------------------------------------------------------
% PURPOSE
%  Compute the element stiffness matrix for two dimensional bar element.
%
% INPUT:  ex = [x1 x2];
%         ey = [y1 y2];      element node coordinates
%
%         ep = [E A]         E: Young's modulus
%                            A: Cross section area
%
% OUTPUT: Ke : stiffness matrix, dim(Ke)= 4 x 4
%----------------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  E=ep(1);  A=ep(2); 

  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);

  Kle=E*A/L*[ 1 -1; 
             -1  1];

  n=b'/L;  G=[   n      zeros(size(n));  
              zeros(size(n))     n   ];

  Ke=G'*Kle*G;
%--------------------------end--------------------------------
