 function [Ke,fe]=beam1e(ex,ep,eq)
% Ke=beam1e(ex,ep)
% [Ke,fe]=beam1e(ex,ep,eq)
%--------------------------------------------------------------------------
%    PURPOSE
%      Compute the stiffness matrix for a one dimensional beam element.
%
%    INPUT:    ex = [x1 x2]    element node coordinates
%              ep = [E I]      element properties;
%                                     E: Young's modulus
%                                     I: moment of inertia
%              eq = [qY]       distributed load 
%            
%    OUTPUT:   Ke : beam stiffness matrix (4 x 4)
%              fe : element load vector (4 x 1)
%--------------------------------------------------------------------------
  
% LAST MODIFIED: O Dahlblom    2019-01-09
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  E=ep(1);  I=ep(2);
  DEI=E*I;

  qY=0;  if nargin>2; qY=eq(1); end

  L=ex(2)-ex(1);

  Ke=DEI/(L^3)*[ 12  6*L   -12  6*L;
                 6*L 4*L^2  -6*L 2*L^2;
                -12 -6*L    12 -6*L;
                 6*L 2*L^2 -6*L 4*L^2];
 
  fe=qY*[L/2 L^2/12 L/2 -L^2/12]';
%------------------------------------end-----------------------------------
  
