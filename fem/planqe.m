 function [Ke,fe]=planqe(ex,ey,ep,D,eq)
% Ke=planqe(ex,ey,ep,D)
% [Ke,fe]=planqe(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a quadrilateral
%  plane stress or plane strain element.
%
% INPUT:    ex=[x1 x2 x3 x4]    element coordinates
%           ey=[y1 y2 y3 y4]
%                             
%           ep = [ptype t]      ptype: analysis type
%                               t: element thickness 
%
%           D                   constitutive matrix
%
%           eq = [bx;           bx: body force in x direction
%                 by]           by: body force in y direction
%
% OUTPUT: Ke :  element stiffness matrix (8 x 8)
%         fe : equivalent nodal forces (8 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-15
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

K=zeros(10,10);   f=zeros(10,1);

xm=sum(ex)/4; ym=sum(ey)/4;

if nargin==4  b1=[0;0]; else  b1=eq; end

[ke1,fe1]=plante([ex(1) ex(2) xm],[ey(1) ey(2) ym],ep,D,b1);
[K,f]=assem([1 1 2 3 4 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante([ex(2) ex(3) xm],[ey(2) ey(3) ym],ep,D,b1);
[K,f]=assem([2 3 4 5 6 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante([ex(3) ex(4) xm],[ey(3) ey(4) ym],ep,D,b1);
[K,f]=assem([3 5 6 7 8 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante([ex(4) ex(1) xm],[ey(4) ey(1) ym],ep,D,b1);
[K,f]=assem([4 7 8 1 2 9 10],K,ke1,f,fe1);
[Ke,fe]=statcon(K,f,[9;10]);

%--------------------------end--------------------------------
