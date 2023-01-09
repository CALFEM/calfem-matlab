 function [Ke,fe]=bar1e(ex,ep,eq)
% Ke=bar1e(ex,ep)
% [Ke,fe]=bar1e(ex,ep,eq)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute the stiffness matrix for a onedimensional bar element.
%
%    INPUT:    ex = [x1 x2]        element node coordinates
% 
%              ep = [E A]   element properties;
%                                     E: Young's modulus
%                                     A: cross section area
%
%              eq = [qx]              distributed load (local direction)
%            
%    OUTPUT:   Ke : bar stiffness matrix (2 x 2)
%
%              fe : element load vector (2 x 1)
%-------------------------------------------------------------------------
  
% LAST MODIFIED: O Dahlblom  2015-10-22
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------

  L=[ex(2)-ex(1)];
%
  E=ep(1);  A=ep(2);
  DEA=E*A
% 
  qx=0;  if nargin>2; qx=eq(1); end
%
  Ke =DEA/L*[1 -1;
            -1  1];
%   
  fe=qx*L*[1/2  1/2]';
%------------------------------------end-------------------------------------
  
