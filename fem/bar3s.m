function [es]=bar3s(ex,ey,ez,ep,ed)
% es=bar3s(ex,ey,ez,ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute normal force in three dimensional bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]       element node coordinates
%         ez = [z1 z2] 
%    
%         ep = [E A]         element properties   
%                               E : Young's modulus
%                               A : Cross section area
%
%         ed : [u1 ... u6]   element displacements
%
% OUTPUT: es = [N]     normal force
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
 b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
 L=sqrt(b'*b);

 n=b'/L;   G=[   n      zeros(size(n));
              zeros(size(n))     n   ];

 E=ep(1); A=ep(2); Kle=E*A/L*[ 1 -1;
                              -1  1];

 N=E*A/L*[-1 1]*G*ed';
 es=N;
%--------------------------end--------------------------------
