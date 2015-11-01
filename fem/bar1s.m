function [es]=bar1s(ep,ed)
% es=bar1s(ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute force in a number of identical (nie) spring elements.
%
% INPUT:  ep = [k]        spring stiffness or analog quantity.
%         ed = [u1 u2;    element displacement vector
%               .....]    u1, u2: nodal displacements
%                            one row for each element
%
% OUTPUT: es  = [N1;      normal force,one row for each element 
%                N2;        
%                .. ]     dim(es)= nie x 1 
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

  k = ep;
  u=ed';
  P1 = (k*[-1 1]*u)';
  es=P1;
%--------------------------end--------------------------------

