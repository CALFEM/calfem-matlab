 function [Ke,fe]=planre(ex,ey,ep,D,eq)
% Ke=planre(ex,ey,ep,D)
% [Ke,fe]=planre(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a rectangular Melosh  
%  plane stress or plane strain element.
%  NOTE! Element sides must be parallell to the coordinate axis.
%
% INPUT:  ex = [x1 x3]           element coordinates
%         ey = [y1 y3]
% 
%         ep = [ptype t ]        ptype: analysis type
%                                t: thickness
%                                    
%         D                      constitutive matrix
%
%         eq = [bx;              bx: body force in x direction
%               by]              by: body force in y direction
%
% OUTPUT: Ke : element stiffness matrix (8 x 8)
%         fe : equivalent nodal forces (8 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa    1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
   
ptype=ep(1);
t=ep(2);

bx=0 ; by=0 ;
if nargin==5; bx=eq(1) ; by=eq(2) ; end

Ke=zeros(8);
a=(ex(2)-ex(1))/2;  b=(ey(2)-ey(1))/2;

xgp=[-1  1 1 -1]/sqrt(3);
ygp=[-1 -1 1  1]/sqrt(3);

%--------- plane stress --------------------------------------
if ptype==1
   
     colD=size(D,2);
     if colD>3
       Cm=inv(D);
       Dm=inv(Cm([1 2 4],[1 2 4]));
     else
       Dm=D;
     end
     
     for i=1:4   
       x=xgp(i)*a; y=ygp(i)*b;
    
       B=[-(b-y)    0     b-y     0   b+y   0  -(b+y)   0  ;
             0   -(a-x)    0   -(a+x)  0   a+x    0    a-x ;
          -(a-x) -(b-y) -(a+x)   b-y  a+x  b+y   a-x -(b+y)]/(4*a*b);
    
       Ke=Ke+B'*Dm*B*a*b*t;
     end

     fe=a*b*[bx by bx by bx by bx by]'*t;

%--------- plane strain --------------------------------------
elseif ptype==2

     colD=size(D,2);
     if colD>3
       Dm=D([1 2 4],[1 2 4]);
     else
       Dm=D;
     end
     
     for i=1:4   
       x=xgp(i)*a; y=ygp(i)*b;
    
       B=[-(b-y)    0     b-y     0   b+y   0  -(b+y)   0  ;
             0   -(a-x)    0   -(a+x)  0   a+x    0    a-x ;
          -(a-x) -(b-y) -(a+x)   b-y  a+x  b+y   a-x -(b+y)]/(4*a*b);
    
       Ke=Ke+B'*Dm*B*a*b*t;
     end

     fe=a*b*[bx by bx by bx by bx by]'*t;

else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
