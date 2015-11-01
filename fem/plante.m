  function [Ke,fe]=plante(ex,ey,ep,D,eq)
% Ke=plante(ex,ey,ep,D)
% [Ke,fe]=plante(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a triangular plane stress
%  or plane strain element.
%
% INPUT:  ex = [x1 x2 x3]         element coordinates
%         ey = [y1 y2 y3]
% 
%         ep = [ptype t ]         ptype: analysis type
%                                 t: thickness
% 
%         D                       constitutive matrix
%
%         eq = [bx;               bx: body force x-dir
%               by]               by: body force y-dir
%
% OUTPUT: Ke : element stiffness matrix (6 x 6)
%         fe : equivalent nodal forces (6 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

ptype=ep(1);
t=ep(2);

bx=0.; by=0.; if nargin==5;  bx=eq(1); by=eq(2); end

C=[ 1  ex(1) ey(1)   0     0       0  
    0    0     0     1   ex(1)   ey(1)
    1  ex(2) ey(2)   0     0       0  
    0    0     0     1   ex(2)   ey(2)
    1  ex(3) ey(3)   0     0       0  
    0    0     0     1   ex(3)   ey(3)];

A=1/2*det([ones(3,1) ex' ey']);

%--------- plane stress --------------------------------------
if ptype==1 
       B=[0 1 0 0 0 0
          0 0 0 0 0 1
          0 0 1 0 1 0]*inv(C);
       
       colD=size(D,2);
       if colD>3
         Cm=inv(D);
         Dm=inv(Cm([1 2 4],[1 2 4]));
       else
         Dm=D;
       end

       Ke=B'*Dm*B*A*t;
       fe=A/3*[bx by bx by bx by]'*t;

%--------- plane strain --------------------------------------       
elseif ptype==2
       B=[0 1 0 0 0 0
          0 0 0 0 0 1
          0 0 1 0 1 0]*inv(C);

       colD=size(D,2);
       if colD>3
         Dm=D([1 2 4],[1 2 4]);
       else
         Dm=D;
       end

       Ke=B'*Dm*B*A*t;
       fe=A/3*[bx by bx by bx by]'*t;
       
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
