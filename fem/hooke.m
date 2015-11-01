  function [D]=hooke(ptype,E,v)
% D=hooke(ptype,E,v)
%-------------------------------------------------------------
%  PURPOSE
%   Calculate the material matrix for a linear
%   elastic and isotropic material.
%
% INPUT:  ptype=1:  plane stress
%               2:  plane strain
%               3:  axisymmetry
%               4:  three dimensional
%
%          E : Young's modulus
%          v : Poissons const.
%
% OUTPUT: D : material matrix
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
 if ptype==1
        Dm=E/(1-v^2)*[1  v   0;
                      v  1   0;
                      0  0 (1-v)/2];
 elseif ptype==2
        Dm=E/(1+v)/(1-2*v)*[1-v  v    v        0;
                             v  1-v   v        0;
                             v   v   1-v       0;
                             0   0    0   (1-2*v)/2];;
 elseif ptype==3
        Dm=E/(1+v)/(1-2*v)*[1-v  v    v        0;
                             v  1-v   v        0;
                             v   v   1-v       0;
                             0   0    0   (1-2*v)/2];;
 elseif ptype==4
        Dm=E/(1+v)/(1-2*v)*[1-v  v    v    0    0    0;
                             v  1-v   v    0    0    0;
                             v   v   1-v   0    0    0;
                             0   0    0 (1-2*v)/2    0    0;
                             0   0    0    0 (1-2*v)/2    0;
                             0   0    0    0    0  (1-2*v)/2];
 else
   error('Error ! Check first argument, ptype=1,2,3 or 4 allowed')
   return
 end
 D=Dm;
%--------------------------end--------------------------------

