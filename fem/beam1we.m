 function [Ke,fe]=beam1we(ex,ep,eq)
% Ke=beam1we(ex,ep)
% [Ke,fe]=beam1we(ex,ep,eq)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute the stiffness matrix for a one dimensional beam element 
%      on elastic foundation.
%
%    INPUT:    ex = [x1 x2]    element node coordinates
%              ep = [E I kY]   element properties;
%                                     E: Young's modulus
%                                     I: moment of inertia
%                                     kY: transversal found. stiffness
%              eq = [qY]       distributed load 
%            
%    OUTPUT:   Ke : beam stiffness matrix (4 x 4)
%              fe : element load vector (4 x 1)
%-------------------------------------------------------------------------
  
% LAST MODIFIED: O Dahlblom    2016-02-17
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
  E=ep(1);  I=ep(2);  kY=ep(3);
  DEI=E*I;

  qY=0;  if nargin>2; qY=eq(1); end

  L=ex(2)-ex(1);
 
  K0=DEI/(L^3)*[ 12  6*L   -12  6*L;
                 6*L 4*L^2  -6*L 2*L^2;
                -12 -6*L    12 -6*L;
                 6*L 2*L^2 -6*L 4*L^2];
         
  Ks=kY*L/420*[156   22*L   54   -13*L ;
            22*L  4*L^2  13*L  -3*L^2;
            54    13*L   156   -22*L ;
            -13*L -3*L^2 -22*L 4*L^2];
   
  Ke=K0+Ks;
 
  fe=qY*[L/2 L^2/12 L/2 -L^2/12]';
%------------------------------------end-------------------------------------
  
