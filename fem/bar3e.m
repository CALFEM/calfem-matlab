function [Ke]=bar3e(ex,ey,ez,ep)
% Ke=bar3e(ex,ey,ez,ep)
%--------------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for three dimensional bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]         element node coordinates
%         ez = [z1 z2]
%    
%         ep = [E A]           E: Young's modulus
%                              A: cross section area
%
% OUTPUT: Ke : stiffness matrix, dim(Ke)= 6 x 6
%--------------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%--------------------------------------------------------------------
  E=ep(1);  A=ep(2); 
 
  b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
  L=sqrt(b'*b);

  Kle=E*A/L*[1 -1;
            -1  1];

  n=b'/L;   G=[   n   zeros(size(n));
                zeros(size(n))   n   ]; 

  Ke=G'*Kle*G;
%--------------------------end--------------------------------
