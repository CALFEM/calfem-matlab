 function [Ke,fe]=plantce(ex,ey,ep,eq)
% Ke=plantce(ex,ey,ep)
% [Ke,fe]=plantce(ex,ey,ep,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a rectangular Turner-Clough
%  plane stress or plane strain element.
%  NOTE! Element sides must be parallel to the coordinate axis.
%
% INPUT:  ex = [x1 x3]           element coordinates
%         ey = [y1 y3]
% 
%         ep = [ptype t E v ]   ptype: 1 -> plane stress
%                                      2 -> plane strain
%                                t: thickness
%                                E: Young's modulus
%                                v: Poisson's ratio
%    
%         eq = [bx;              bx: body force in x direction
%               by]              by: body force in y direction
%
% OUTPUT: Ke : element stiffness matrix (8 x 8)
%         fe : equivalent nodal forces (8 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: A Olsson    2002-12-16
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

ptype=ep(1) ; t=ep(2); E=ep(3); v=ep(4);
 
bx=0 ; by=0 ;   
if nargin==4; bx=eq(1) ; by=eq(2) ; end
   
Ke=zeros(8);
a=(ex(2)-ex(1))/2  ;  b=(ey(2)-ey(1))/2  ;

xgp=[-1 1 1 -1]/sqrt(3);
ygp=[-1 -1 1 1]/sqrt(3); 
 
%--------- plane stress --------------------------------------
if ptype==1
   
   D=hooke(ptype,E,v);

   for i=1:4
     x=xgp(i)*a; y=ygp(i)*b;
     B=[-(b-y) -v*x   b-y   v*x   b+y  -v*x -(b+y)  v*x ;
         -v*y -(a-x)  v*y -(a+x) -v*y   a+x   v*y   a-x ;
          -a    -b    -a     b     a     b     a    -b  ]/(4*a*b) ;
     Ke=Ke+B'*D*B*a*b*t;
   end

   fe=a*b*[bx by bx by bx by bx by]'*t;
 
%--------- plane strain --------------------------------------
elseif ptype==2

   D=hooke(ptype,E,v);
   D=D([1 2 4],[1 2 4]);
   
   for i=1:4
     x=xgp(i)*a; y=ygp(i)*b;
     B=[-(b-y) -v*x   b-y   v*x   b+y  -v*x -(b+y)  v*x ;
         -v*y -(a-x)  v*y -(a+x) -v*y   a+x   v*y   a-x ;
          -a    -b    -a     b     a     b     a    -b  ]/(4*a*b) ;
     Ke=Ke+B'*D*B*a*b*t;
   end

   fe=a*b*[bx by bx by bx by bx by]'*t;
 
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
