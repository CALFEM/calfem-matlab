function [he]=cp2s3f(ex,ey,ep)
% [he]=cp2s3f(ex,ey,ep)
%-----------------------------------------------------------
% PURPOSE
%  Compute element coupling matrix between a 8 node 
%  isoparametric acoustic element and a 2 node beam element. 
%  A non-symmetric pressure
%  formulation according to Eq. (3.46) in Carlsson TVSM-1005
%  is used. Note that the corresponding coupling matrices 
%  used in the non-symmetric Psi-formulation according to 
%  Eq. (3.47) is just the transpose of the matrices below. 
%
%     1-------------2 Beam    
%     *-------------*
%     1------2------3 Fluid
%     |             |
%     |             |
% 
% INPUT:  ex = [x1 x2]   element coordinates
%         ey = [y1 y2]
%                             
%         ep = [t]       thickness 
%
% OUTPUT: he :  element coupling matrix (6 x 3)
%-----------------------------------------------------------

% LAST MODIFIED: G Sandberg    1996-03-08
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

t=ep(1); 

b=[ex(2)-ex(1);
   ey(2)-ey(1)];
   
L=sqrt(b'*b);

detJ=0.5*L;
if detJ<10*eps
  disp('Jacobideterminanten lika med noll!')
end

% Loop over gauss points (ir=3)

Cle=zeros(6,3);

for ip=1:3
  if ip==1
    xa= -0.774596699241483;
    xh=  0.555555555555555;
  elseif ip==2
    xa= 0.;
    xh= 0.888888888888888;
  elseif ip==3
    xa=  0.774596699241483;
    xh=  0.555555555555555;
  end 

  % Structural shape functions

  a=xa+1;
  s(1)=0;
  s(2)=(4-3*a*a+a*a*a)*0.25;
  s(3)=(a-a*a+a*a*a*0.25)*0.5*L;
  s(4)=0;
  s(5)=(3*a*a-a*a*a)*0.25;
  s(6)=-(a*a-a*a*a*0.5)*0.25*L;

  % Fluid shape functions

  p(1)=0.5*xa*(xa-1);
  p(2)=1-xa*xa;
  p(3)=0.5*xa*(xa+1);

  for i=1:6; for j=1:3
      Cle(i,j)=Cle(i,j)+s(i)*p(j)*xh*detJ;
  end; end
 
end

% Local to global transformation matrix G 

n=b/L;

G=[n(1) n(2) 0    0    0   0;
  -n(2) n(1) 0    0    0   0;
   0    0    1    0    0   0;
   0    0    0   n(1) n(2) 0;
   0    0    0  -n(2) n(1) 0;
   0    0    0    0    0   1];

Ce=G'*Cle;

he= t*Ce;
%-------------------------- end -------------------------------
