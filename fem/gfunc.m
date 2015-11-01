function [t,g]=gfunc(G,dt);
% [t,g]=gfunc(G,dt)
%-------------------------------------------------------------
% PURPOSE
%  Form vector with function values at equally spaced  
%  points by linear interpolation. 
%
% INPUT:  G = [t_i g_i]     t_i: time i
%                           g_i: g(t_i)
%                           dim(G)=np x 2, np= number of points
%  
% OUTPUT: t : 1-D vector with equally spaced time points 
%         g : 1-D vector with corresponding function values 
%-------------------------------------------------------------

% LAST MODIFIED: K Persson 1997-04-07
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [np nc]=size(G);
  ti=G(1,1):dt:G(np,1);
  g1=interp1(G(:,1),G(:,2),ti');
%-------------------------------------------------------------
if nargout == 0
  disp('Discrete function values have been computed');
  disp(g1);
  return
end
t=ti;
g=g1';
%--------------------------end--------------------------------
