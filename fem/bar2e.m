function [Ke,fe]=bar2e(ex,ey,ep,eq)
% Ke=bar2e(ex,ey,ep)
% [Ke,fe]=bar2e(ex,ey,ep,eq)
%----------------------------------------------------------------------
% PURPOSE
%  Compute the element stiffness matrix for two dimensional bar element.
%
% INPUT:  ex = [x1 x2];
%         ey = [y1 y2];      element node coordinates
%
%         ep = [E A]         E: Young's modulus
%                            A: Cross section area
%         eq = [qX]          distributed load, local direction
%
% OUTPUT: Ke : stiffness matrix (4 x 4)
%         fe : element load vector (4 x 1)
%----------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2015-10-20
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  E=ep(1);  A=ep(2); 
  DEA=E*A;
  
  qX=0; if nargin>3; qX=eq(1); end

  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  L=sqrt(dx*dx+dy*dy);

  Kle=DEA/L*[ 1 -1; 
             -1  1];

  fle=L*[qX/2 qX/2]';

  nxX=dx/L;
  nyX=dy/L;
  G=[nxX nyX  0   0 ;  
      0   0  nxX nyX];

  Ke=G'*Kle*G;   fe=G'*fle;
%--------------------------end--------------------------------
