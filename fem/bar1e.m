function [Ke]=bar1e(ep);
% Ke=bar1e(ep)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix
%  for spring (analog) element.
%
% INPUT:  ep = [k];       spring stiffness or analog quantity.
%
% OUTPUT: Ke :            stiffness matrix, dim(Ke)= 2 x 2
%-------------------------------------------------------------

% LAST MODIFIED: K Persson   1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
k = ep;
Ke1 = [ k -k;
       -k  k];
Ke=Ke1;
%--------------------------end--------------------------------
