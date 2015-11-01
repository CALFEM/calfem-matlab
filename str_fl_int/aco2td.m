function [Ke,Me,fe]=aco2td(ex,ey,ep,eq)
% [Ke,Me]=aco2td(ex,ey,ep)
% [Ke,Me,fe]=aco2td(ex,ey,ep,eq)
%----------------------------------------------------------
% PURPOSE
%  Compute element stiffness and consistent element
%  mass matrices for the triangular acoustic element.
%
% INPUT:  ex = [x1 x2 x3]
%         ey = [y1 y2 y3]      element coordinates
%
%         ep = [t c raa]       thickness,speed of sound and 
%                              density      
%                             
%         eq                  mass inflow per unit volume and time
%                              (second derivative)  
%
% OUTPUT: Ke :       element stiffness matrix (3 x 3)
%         Me :       element mass matrix (3 x 3)
%         fe :       element load vector (3 x 1)
%----------------------------------------------------------

% LAST MODIFIED: G Sandberg    1996-03-09
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%----------------------------------------------------------
  t=ep(1); c=ep(2); raa=ep(3);
  if nargin==3; qe=0; end 
%
  C=[ones(3,1) ex' ey'];   B=[0 1 0;
                              0 0 1 ]*inv(C);   A=1/2*det(C);
  Ke1=t*c*c*B'*B*A;
%
  Me1=t*A/12*[2 1 1;
              1 2 1;
              1 1 2];
%                   
  fe1=t*c*c*qe*A/3*[1 1 1]';

  Ke=Ke1;  Me=Me1;  fe=fe1; 
%------------------------- end -----------------------------
