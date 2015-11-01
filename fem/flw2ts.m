function [es,et]=flw2ts(ex,ey,D,ed)
% [es,et]=flw2ts(ex,ey,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute flows or corresponding quantities in the
%  triangular field element.
%
% INPUT:  ex = [x1 x2 x3]
%         ey = [y1 y2 y3]         element coordinates
%                             
%         D = [kxx kxy;
%              kyx kyy]           constitutive matrix
%
%         ed =[u1 u2 u3]          u1,u2,u3: nodal values
%              .. .. ..;
%
%
% OUTPUT: es=[ qx qy ] 
%              ... ..]            element flows
%
%         et=[ gx gy ]
%              ... ..]            element gradients
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1997-04-14
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  C = [ones(3,1) ex' ey' ];
  
  B = [0 1 0;
       0 0 1 ]*inv(C);

  qs=-D*B*ed';
  qt=B*ed';

  es=[qs'];
  et=[qt'];
%--------------------------end--------------------------------
