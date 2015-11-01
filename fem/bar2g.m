function [Ke]=bar2g(ex,ey,ep,N)
% Ke=bar2g(ex,ey,ep,N)
%----------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for two dimensional geometric
%  nonlinear bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]       element node coordinates
%
%         ep = [E A]         E: Young's modulus
%                            A: Cross section area
%
%         N                  normal force
%
% OUTPUT: Ke :   stiffness matrix, dim(Ke)= 4 x 4
%----------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%----------------------------------------------------------------
   E=ep(1);   A=ep(2); 

   b=[ ex(2)-ex(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);

   n=b'/L;  G=[n(1) n(2)  0   0   ;
              -n(2) n(1)  0   0   ;
                0    0   n(1) n(2);
                0    0  -n(2) n(1)];

   Kle=E*A/L*[1  0 -1  0;
              0  0  0  0;
             -1  0  1  0;
              0  0  0  0] + N/L*[0  0  0  0;
                                 0  1  0 -1;
                                 0  0  0  0;
                                 0 -1  0  1];

   Ke=G'*Kle*G;
%--------------------------end--------------------------------

