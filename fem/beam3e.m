 function [Ke,fe]=beam3e(ex,ey,ez,eo,ep,eq)
% Ke=beam3e(ex,ey,ez,eo,ep)
% [Ke,fe]=beam3e(ex,ey,ez,eo,ep,eq)
%----------------------------------------------------------------
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
%            eq = [qx qy qz qw];    distributed loads
%
%    OUTPUT: Ke : beam stiffness matrix (12 x 12)
%
%            fe : equivalent nodal forces (12 x 1)
%-----------------------------------------------------------------  

% LAST MODIFIED: E Serrano    1995-09-21 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------


  b=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
  L=sqrt(b'*b);  n1=b/L;

  lc=sqrt(eo*eo'); n3=eo/lc;
    
 %
    
     if nargin==5;   eq=[0 0 0 0];  end 
   
     qx=eq(1); qy=eq(2); qz=eq(3); qw=eq(4);
  %
    E=ep(1); Gs=ep(2);
    A=ep(3);
    Iy=ep(4); Iz=ep(5);
    Kv=ep(6);
 
    a=E*A/L       ; b=12*E*Iz/L^3 ; c=6*E*Iz/L^2;
    d=12*E*Iy/L^3 ; e=6*E*Iy/L^2  ; f=Gs*Kv/L;
    g=2*E*Iy/L    ; h=2*E*Iz/L    ;

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
   fle=L/2*[qx qy qz qw -1/6*qz*L 1/6*qy*L qx qy qz qw 1/6*qz*L -1/6*qy*L]';

 %
    n2(1)=n3(2)*n1(3)-n3(3)*n1(2);
    n2(2)=-n1(3)*n3(1)+n1(1)*n3(3);
    n2(3)=n3(1)*n1(2)-n1(1)*n3(2);
%
    An=[n1';
        n2;
        n3];
%
    G=[  An     zeros(3) zeros(3) zeros(3);
       zeros(3)   An     zeros(3) zeros(3);
       zeros(3) zeros(3)   An     zeros(3);
       zeros(3) zeros(3) zeros(3)   An    ];
%


 %
    Ke1=G'*Kle*G;  fe1=G'*fle;
 %----------------------------------------------------------   
    if nargout==0
      disp('Element stiffness matrix: ');
      disp(Ke1);
      
      if nargin==6
        disp('Element load vector: ');
        disp(fe1);
      end
      return 
    end
    Ke=Ke1;
  
    if nargin==6  fe=fe1;  end
 %-------------------------- end -------------------------------

