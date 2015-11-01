function [es]=bar2s(ex,ey,ep,ed)
% es=bar2s(ex,ey,ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute normal force in two dimensional bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]        element coordinates
%
%         ep = [E A]          E : Young's modulus
%                             A : Cross section area
%
%         ed : [u1 u2 u3 u4]  element displacement vector
%
% OUTPUT: es = [N]            element force 
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
   E=ep(1);  A=ep(2);
 
   b=[ ex(2)-ex(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);

   Kle=E*A/L*[1 -1 ;
             -1  1 ];

   n=b'/L;   G=[   n      zeros(size(n));
                zeros(size(n))     n    ];

   u=ed';
   N=E*A/L*[-1 1]*G*u;
   es=N;
%--------------------------end--------------------------------
