 function [Ke,fe]=beam3e(ex,ey,ez,eo,ep,eq)
% Ke=beam3e(ex,ey,ez,eo,ep)
% [Ke,fe]=beam3e(ex,ey,ez,eo,ep,eq)
%--------------------------------------------------------------------------
%    PURPOSE
%       Calculate the stiffness matrix for a 3D elastic Bernoulli
%       beam element. 
% 
%    INPUT:  ex = [x1 x2]        
%            ey = [y1 y2]   
%            ez = [z1 z2]           node coordinates  
%
%            eo = [xz yz zz];       orientation of local z axis
%
%            ep = [E G A Iy Iz Kv]; element properties 
%                                   E: Young's modulus
%                                   G: Shear modulus 
%                                   A: Cross section area
%                                   Iy: moment of inertia,local y-axis
%                                   Iz: moment of inertia,local z-axis
%                                   Kv: Saint-Venant's torsion constant
% 
%            eq = [qX qY qZ qW];    distributed loads
%
%    OUTPUT: Ke : beam stiffness matrix (12 x 12)
%
%            fe : equivalent nodal forces (12 x 1)
%--------------------------------------------------------------------------  

% LAST MODIFIED: O Dahlblom    2015-10-19 
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------

  E=ep(1); Gs=ep(2);
  A=ep(3); Iy=ep(4); Iz=ep(5); Kv=ep(6);
  DEA=E*A; DEIz=E*Iz; DEIy=E*Iy; DGK=Gs*Kv;

  qX=0; qY=0; qZ=0; qW=0; 
  if nargin>5; qX=eq(1); qY=eq(2); qZ=eq(3); qW=eq(4); end
  
  dx=ex(2)-ex(1);
  dy=ey(2)-ey(1);
  dz=ez(2)-ez(1);
  L=sqrt(dx*dx+dy*dy+dz*dz);
  n1=[dx dy dz]/L;
  lc=sqrt(eo*eo'); 
  n3=eo/lc;
  
  if nargin==5;   eq=[0 0 0 0];  end 
   
  a=DEA/L       ; b=12*DEIz/L^3 ; c=6*DEIz/L^2;
  d=12*DEIy/L^3 ; e=6*DEIy/L^2  ; f=DGK/L;
  g=2*DEIy/L    ; h=2*DEIz/L    ;

  Kle=[ a  0  0  0  0  0 -a  0  0  0  0  0 ;
        0  b  0  0  0  c  0 -b  0  0  0  c ;
        0  0  d  0 -e  0  0  0 -d  0 -e  0 ;
        0  0  0  f  0  0  0  0  0 -f  0  0 ;
        0  0 -e  0 2*g 0  0  0  e  0  g  0 ;
        0  c  0  0  0 2*h 0 -c  0  0  0  h ;
       -a  0  0  0  0  0  a  0  0  0  0  0 ;
        0 -b  0  0  0 -c  0  b  0  0  0 -c ;
        0  0 -d  0  e  0  0  0  d  0  e  0 ;
        0  0  0 -f  0  0  0  0  0  f  0  0 ;
        0  0 -e  0  g  0  0  0  e  0 2*g 0 ;
        0  c  0  0  0  h  0 -c  0  0  0 2*h];
 %
  fle=L*[qX/2 qY/2 qZ/2 qW/2 -qZ*L/12 qY*L/12 qX/2 qY/2 qZ/2 qW/2 qZ*L/12 -qY*L/12]';

 %
  n2(1)= n3(2)*n1(3)-n3(3)*n1(2);
  n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
  n2(3)= n3(1)*n1(2)-n1(1)*n3(2);
%
  An=[n1; n2; n3];
%
  G=[  An     zeros(3) zeros(3) zeros(3);
     zeros(3)   An     zeros(3) zeros(3);
     zeros(3) zeros(3)   An     zeros(3);
     zeros(3) zeros(3) zeros(3)   An    ];
%
    Ke=G'*Kle*G;  fe=G'*fle;
%--------------------------- end ------------------------------------------

