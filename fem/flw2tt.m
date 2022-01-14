function [Ke,Ce,fe]=flw2tt(ex,ey,ep,D,eq)
% Ke=flw2te(ex,ey,ep,D)
% [Ke,fe]=flw2te(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness (conductivity) and capacitivity 
%  matrices for a transient triangular field element.
%
% INPUT:  ex = [x1 x2 x3]
%         ey = [y1 y2 y3]      element coordinates
%
%         ep = [t, rho, c]     element thickness, density, capacitivity  	 
%                             
%         D = [kxx kxy;
%              kyx kyy]        constitutive matrix
%
%         eq                   heat supply per unit volume
%
% OUTPUT: Ke :  element 'stiffness' matrix (3 x 3)
%
%         fe :  element load vector (3 x 1)
%-------------------------------------------------------------
% LAST MODIFIED: K Persson   1995-08-23
% Copyright (c) 1991-94 by Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  t=ep(1); rho=ep(2); c=ep(3);
  if nargin==4; eq=0; end 
%
  C=[ones(3,1) ex' ey'];   B=[0 1 0;
                              0 0 1 ]*inv(C);   A=1/2*det(C);
  Ce1=[1/6 1/12 1/12;
       1/12 1/6 1/12;
       1/12 1/12 1/6];
   
  Ke1=B'*D*B*t*A;
  
  fe1=eq*A*t/3*[1 1 1]'; 

  Ce=rho*c*A*t*Ce1;
  Ke=Ke1;  fe=fe1;
%--------------------------end--------------------------------
