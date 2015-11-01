function [Ke]=spring1e(ep);
% Ke=spring1e(ep)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness matrix for spring element.
%
% INPUT:  ep = [k];       spring stiffness or analog quantity
%
% OUTPUT: Ke :            stiffness matrix, dim(Ke)= 2 x 2
%-------------------------------------------------------------

% LAST MODIFIED: P-E Austrell 1994-11-02
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
k = ep;  
Ke = [ k -k;
      -k  k];
%--------------------------end--------------------------------
