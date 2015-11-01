  function [es,et]=plantcs(ex,ey,ep,ed)
% [es,et]=plancs(ex,ey,ep,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a rectangular 
%  Turner-Clough plane stress or plane strain element.
%  NOTE! Element sides must be parallell to the coordinate axis.
%
% INPUT:  ex = [x1 x3]           element coordinates
%         ey = [y1 y3]
% 
%         ep = [ ptype t E v ]   ptype: 1 -> plane stress
%                                       2 -> plane strain
%                                t: thickness
%                                E: Young's modulus
%                                v: Poisson's ratio
%    
%         ed = [u1 u2 .. u8;        element displacement vector
%              ........... ]        one row for each element
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy    element stress matrix
%                ......                ]   one row for each element
%         et = [ epsx epsy [epsz] gamxy    element strain matrix
%                  ......              ]   one row for each element
%-------------------------------------------------------------

% LAST MODIFIED: A Olsson    2002-12-16
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

ptype=ep(1); t=ep(2); E=ep(3); v=ep(4); 

rowed=size(ed);
rowex=size(ex);  

%--------- plane stress --------------------------------------
if ptype==1
   
   D=hooke(2,E,v);        
   Cm=inv(D);
   Dm=inv(Cm([1 2 4],[1 2 4]));

   if rowex==1 incie=1; else incie=0; end

   x=0;  y=0;
   et=[]; es=[]; ie=1;
   for i=1:rowed
     a=(ex(ie,2)-ex(ie,1))/2;  b=(ey(ie,2)-ey(ie,1))/2;

     B=[-(b-y) -v*x   b-y   v*x   b+y  -v*x -(b+y)  v*x ;
         -v*y -(a-x)  v*y -(a+x) -v*y   a+x   v*y   a-x ;
          -a    -b    -a     b     a     b     a    -b  ]/(4*a*b) ;
  
      ee=B*ed(i,:)';
      ss([1 2 4],1)=Dm*ee;
      ee=Cm*ss;
   
      et=[et; ee'];
      es=[es; ss'];
      
      ie=ie+incie;
   end
   
%--------- plane strain --------------------------------------
elseif ptype==2

   D=hooke(ptype,E,v);        
 
   if rowex==1 incie=1; else incie=0; end

   x=0;  y=0;
   et=[]; es=[];ie=1;
   for i=1:rowed
     a=(ex(ie,2)-ex(ie,1))/2;  b=(ey(ie,2)-ey(ie,1))/2;

     B=[-(b-y) -v*x   b-y   v*x   b+y  -v*x -(b+y)  v*x ;
         -v*y -(a-x)  v*y -(a+x) -v*y   a+x   v*y   a-x ;
          -a    -b    -a     b     a     b     a    -b  ]/(4*a*b) ;
      
     e=B*ed(i,:)';
     ee([1 2 4],1)=e;
     ss=D*ee;
   
     et=[et; ee'];
     es=[es; ss'];
     
     ie=ie+incie;
   end

else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
