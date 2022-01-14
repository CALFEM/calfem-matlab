function [Ke]=bar2ge(ex,ey,ep,QX)
% Ke=bar2ge(ex,ey,ep,QX)
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
%         QX:                axial force in the bar
%
% OUTPUT: Ke :   stiffness matrix, dim(Ke)= 4 x 4
%----------------------------------------------------------------

% LAST MODIFIED: O Dahlblom    2015-12-17
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%----------------------------------------------------------------
   E=ep(1);   A=ep(2); 
   DEA=E*A;

   dx=ex(2)-ex(1);
   dy=ey(2)-ey(1);
   L=sqrt(dx*dx+dy*dy);

   K0le=DEA/L*[1  0 -1  0;
               0  0  0  0;
              -1  0  1  0;
               0  0  0  0]; 
   Ksle=QX/L*[0  0  0  0;
              0  1  0 -1;
              0  0  0  0;
              0 -1  0  1];
   Kle=K0le+Ksle;
   
   nxX=dx/L;
   nyX=dy/L;
   nxY=-dy/L;
   nyY=dx/L;
   G=[nxX nyX  0   0 ;  
      nxY nyY  0   0;
       0   0  nxX nyX
       0   0  nxY nyY];
 
   Ke=G'*Kle*G;
%--------------------------end--------------------------------

