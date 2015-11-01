 function [es,et]=planqs(ex,ey,ep,D,ed,eq)
% [es,et]=planqs(ex,ey,ep,D,ed,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a quadrilateral 
%  plane stress or plane strain element.
%
% INPUT:  ex = [x1 x2 x3 x4]      element coordinates
%         ey = [y1 y2 y3 y4]
%
%         ep = [ptype t ]         ptype: analysis type
%                                 t:  thickness
%                                
%         D                       constitutive matrix
%
%         ed = [u1 u2 ..u8;       element displacement vector
%               ..........]       one row for each element
%
%         eq = [bx;               bx: body force in x direction
%               by]               by: body force in y direction
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy    element stress matrix
%                  ......              ]   one row for each element
%         et = [ epsx epsy [epsz] gamxy    element strain matrix
%                  ......              ]   one row for each element
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa   1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

if size(ex,1)~=1 
  error('Error ! PLANQS: only one element at the time');
  return
end

K=zeros(10,10);   f=zeros(10,1);

xm=sum(ex)/4; ym=sum(ey)/4;

if nargin==5  b1=[0;0]; else  b1=eq; end

ex1=[ex(1) ex(2) xm]; ey1=[ey(1) ey(2) ym];
ex2=[ex(2) ex(3) xm]; ey2=[ey(2) ey(3) ym];
ex3=[ex(3) ex(4) xm]; ey3=[ey(3) ey(4) ym];
ex4=[ex(4) ex(1) xm]; ey4=[ey(4) ey(1) ym];

[ke1,fe1]=plante(ex1,ey1,ep,D,b1);
[K,f]=assem([1 1 2 3 4 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante(ex2,ey2,ep,D,b1);
[K,f]=assem([1 3 4 5 6 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante(ex3,ey3,ep,D,b1);
[K,f]=assem([1 5 6 7 8 9 10],K,ke1,f,fe1);
[ke1,fe1]=plante(ex4,ey4,ep,D,b1);
[K,f]=assem([1 7 8 1 2 9 10],K,ke1,f,fe1);

A1=1/2*det([ones(3,1) ex1' ey1']);
A2=1/2*det([ones(3,1) ex2' ey2']);
A3=1/2*det([ones(3,1) ex3' ey3']);
A4=1/2*det([ones(3,1) ex4' ey4']);
Atot=A1+A2+A3+A4;

ni=size(ed,1);
a=[];
for i=1:ni
  a=[a solveq(K,f,[ [1:8]' ed(i,:)'])];
end

[s1,t1]=plants(ex1,ey1,ep,D,[a([1 2 3 4 9 10],:)']);
[s2,t2]=plants(ex2,ey2,ep,D,[a([3 4 5 6 9 10],:)']);
[s3,t3]=plants(ex3,ey3,ep,D,[a([5 6 7 8 9 10],:)']);
[s4,t4]=plants(ex4,ey4,ep,D,[a([7 8 1 2 9 10],:)']);

es=(s1*A1+s2*A2+s3*A3+s4*A4)/Atot;
et=(t1*A1+t2*A2+t3*A3+t4*A4)/Atot;

%--------------------------end--------------------------------
