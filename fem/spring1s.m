function [es]=spring1s(ep,ed)
% es=spring1s(ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute element force in spring element (spring1e).
%
% INPUT:  ep = [k]        spring stiffness or analog quantity
%         ed = [u1 u2]    element displacements
%                         u1, u2: nodal displacements
%                            
%
% OUTPUT: es  = [N]       element force
%-------------------------------------------------------------

% LAST MODIFIED: P-E AUSTRELL 1994-11-02
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  k = ep;  
  u=ed;
  es = k*(u(2)-u(1));
%--------------------------end--------------------------------

