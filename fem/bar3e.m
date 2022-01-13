function [Ke,fe]=bar3e(ex,ey,ez,ep,eq)
% Ke=bar3e(ex,ey,ez,ep)
% [Ke,fe]=bar3e(ex,ey,ez,ep,eq)
%--------------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for three dimensional bar element.
%
% INPUT:  ex = [x1 x2]
%         ey = [y1 y2]        
%         ez = [z1 z2]       element node coordinates
%         ep = [E A]         E: Young's modulus
%                            A: cross section area
%         eq = [qx]          distributed load, local direction
%
% OUTPUT: Ke : stiffness matrix, dim(Ke)= 6 x 6
%         fe : element load vector (6 x 1)
%--------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2015-10-19
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------
  E=ep(1);  A=ep(2); 
  DEA=E*A 

  qX=0; if nargin>4; qX=eq(1); end
 
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  dz=ez(2)-ez(1);
  L=sqrt(dx*dx+dy*dy+dz*dz);

  Kle=DEA/L*[1 -1;
            -1  1];
 
        fle=L*[qX/2 qX/2]';

  nxX=dx/L;
  nyX=dy/L;
  nzX=dz/L;
  G=[nxX nyX nzX  0   0  0 ;  
      0   0   0 nxX nyX nzX];

  Ke=G'*Kle*G;   fe=G'*fle;
%--------------------------end--------------------------------
