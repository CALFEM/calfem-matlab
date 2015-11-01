function [Ke,fe]=flw2qe(ex,ey,ep,D,eq)
% Ke=flw2qe(ex,ey,ep,D)
% [Ke,fe]=flw2qe(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness (conductivity) matrix for a
%  quadrilateral field element, composed of four triangular
%  elements.
%
% INPUT:  ex = [x1 x2 x3 x4]
%         ey = [y1 y2 y3 y4]   element coordinates
%
%         ep = [t]             element thickness                             
%
%         D = [kxx kxy;
%              kyx kyy]        constitutive matrix
%
%         eq                   heat supply per unit volume
%
% OUTPUT: Ke :  element 'stiffness' matrix (4 x 4)
%         fe :  element load vector (4 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
xc=sum(ex)/4;  yc=sum(ey)/4;

K=zeros(5);  f=zeros(5,1);
if nargin==4, q=0; else, q=eq; end

  [k1,f1]=flw2te([ex(1) ex(2) xc],[ey(1) ey(2) yc],ep,D,q);
  [K,f]=assem([1 1 2 5],K,k1,f,f1);
  [k1,f1]=flw2te([ex(2) ex(3) xc],[ey(2) ey(3) yc],ep,D,q);
  [K,f]=assem([2 2 3 5],K,k1,f,f1);
  [k1,f1]=flw2te([ex(3) ex(4) xc],[ey(3) ey(4) yc],ep,D,q);
  [K,f]=assem([3 3 4 5],K,k1,f,f1); 
  [k1,f1]=flw2te([ex(4) ex(1) xc],[ey(4) ey(1) yc],ep,D,q);
  [K,f]=assem([4 4 1 5],K,k1,f,f1);
  [Ke1,fe1]=statcon(K,f,[5]);

  Ke=Ke1; fe=fe1;
%--------------------------end--------------------------------
