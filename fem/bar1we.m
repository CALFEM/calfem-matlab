 function [Ke,fe]=bar1we(ex,ep,eq)
% Ke=bar1we(ex,ep)
% [Ke,fe]=bar1we(ex,ep,eq)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute the stiffness matrix for a onedimensional bar element 
%      with elastic support.
%
%    INPUT:    ex = [x1 x2]        element node coordinates
% 
%              ep = [E A ka]   element properties;
%                                     E: Young's modulus
%                                     A: cross section area
%                                     ka: axial foundation stiffness
%
%              eq = [qx]              distributed load (local direction)
%            
%    OUTPUT:   Ke : beam stiffness matrix (4 x 4)
%
%              fe : element load vector (4 x 1)
%-------------------------------------------------------------------------
  
% LAST MODIFIED: O Dahlblom    2015-12-17
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------

  L=[ex(2)-ex(1)];
%
  E=ep(1);  A=ep(2);  ka=ep(3);
  DEA=E*A;
  % 
  qx=0;  if nargin>2; qx=eq(1); end
%
  K1 =DEA/L*[ 1  -1;
             -1   1];
%         
  K2=ka*L/6*[2   1;
             1   2];
%   
  Ke=K1+K2;
%   
  fe=qx*L*[1/2  1/2]';
%------------------------------------end-------------------------------------
  
