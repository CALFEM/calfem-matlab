  function [es,et]=planrs(ex,ey,ep,D,ed)
% [es,et]=planrs(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for the rectangular 
%  plane stress or plane strain element.
%  NOTE! Element sides must be parallell to the coordinate axis.
%
% INPUT: ex = [x1 x3]                element coordinates
%        ey = [y1 y3]
% 
%        ep = [ptype t ]             ptype: analysis type
%                                    t: thickness
%                                    
%        D                           constitutive matrix
%
%        ed = [u1 u2 .. u8;          element displacement vector
%              ...........]          one row for each element
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy    element stress matrix
%                ......                ]   one row for each element
%         et = [ epsx epsy [gamz] gamxy    element strain matrix
%                  ......              ]   one row for each element
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa   1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
ptype=ep(1);
t=ep(2);
rowed=size(ed,1);
rowex=size(ex,1);

%---------plane stress ---------------------------------------
if ptype==1

   colD=size(D,2);
   if colD>3
     Cm=inv(D);
     Dm=inv(Cm([1 2 4],[1 2 4]));
   else
     Dm=D;
   end
   
   if rowex==1 incie=0; else incie=1; end
   
   x=0 ;  y=0  ;
   et=[];es=[];ie=1;
   for i=1:rowed
   
     a=(ex(ie,2)-ex(ie,1))/2  ;  b=(ey(ie,2)-ey(ie,1))/2  ;

     B=[-(b-y)    0     b-y     0   b+y   0  -(b+y)   0  ;
           0   -(a-x)    0   -(a+x)  0   a+x    0    a-x ;
        -(a-x) -(b-y) -(a+x)   b-y  a+x  b+y   a-x -(b+y)]/(4*a*b);

     ee=B*ed(i,:)';
     if colD>3
       ss=zeros(colD,1);
       ss([1 2 4])=Dm*ee;
       ee=Cm*ss;
     else
       ss=Dm*ee;
     end
     
     et=[et; ee'];
     es=[es; (D*ee)'];
     
     ie=ie+incie;
     
   end
   
%--------- plane strain --------------------------------------
elseif ptype==2

   colD=size(D,2);
   if rowex==1 incie=0; else incie=1; end

   x=0 ;  y=0  ;
   et=[];es=[];ie=1;ee=zeros(colD,1);
   for i=1:rowed
   
     a=(ex(i,2)-ex(i,1))/2  ;  b=(ey(i,2)-ey(i,1))/2  ;

     B=[-(b-y)    0     b-y     0   b+y   0  -(b+y)   0  ;
           0   -(a-x)    0   -(a+x)  0   a+x    0    a-x ;
        -(a-x) -(b-y) -(a+x)   b-y  a+x  b+y   a-x -(b+y)]/(4*a*b);
    
     e=B*ed(i,:)';
     if colD>3 ee([1 2 4])=e; else ee=e; end
     
     et=[et; ee'];
     es=[es; (D*ee)'];
   end

else
  error('Error ! Check first argument, ptype=1 or 2 allowed')
  return
end
%--------------------------end--------------------------------
