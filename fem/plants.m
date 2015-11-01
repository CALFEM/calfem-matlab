  function [es,et]=plants(ex,ey,ep,D,ed)
% [es,et]=plants(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a
%  triangular plane stress or plane strain element.
%
% INPUT:  ex = [x1 x2 x3]         element coordinates
%         ey = [y1 y2 y3]
%
%         ep = [ptype t ]         ptype: analysis type
%                                 t: thickness
% 
%         D                       constitutive matrix
%
%         ed =[u1 u2 ...u6        element displacement vector
%              ......     ]       one row for each element
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy   element stress matrix
%               ......                 ]  one row for each element
%
%         et = [ epsx epsy [epsz] gamxy   element strain matrix
%               ......                 ]  one row for each element
%-------------------------------------------------------------

% LAST MODIFIED: A Olsson 1999-03-01
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
ptype=ep(1);

rowed=size(ed,1);
rowex=size(ex,1);

%--------- plane stress --------------------------------------
if ptype==1 
  
  colD=size(D,2);
  if colD>3 
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end
  
  if rowex==1 incie=0; else incie=1; end
     
  et=[];es=[];ie=1;
  for i=1:rowed
    
    C=[ 1  ex(ie,1) ey(ie,1)   0          0          0  
        0         0        0   1   ex(ie,1)   ey(ie,1)
        1  ex(ie,2) ey(ie,2)   0          0          0  
        0         0        0   1   ex(ie,2)   ey(ie,2)
        1  ex(ie,3) ey(ie,3)   0          0          0  
        0         0        0   1   ex(ie,3)   ey(ie,3)];

    B=[0 1 0 0 0 0;
       0 0 0 0 0 1;
       0 0 1 0 1 0]*inv(C);
    
    ee=B*ed(i,:)';
    if colD>3
      ss=zeros(colD,1);
      ss([1 2 4])=Dm*ee;
      ee=Cm*ss;
    else
      ss=Dm*ee;
    end

    et=[et;ee'];
    es=[es;ss'];

    ie=ie+incie;
  end

%--------- plane strain --------------------------------------
elseif ptype==2

  colD=size(D,2);
  if rowex==1 incie=0; else incie=1; end

  et=[];es=[];ie=1;ee=zeros(colD,1);
  for i=1:rowed
    
    C=[ 1  ex(ie,1) ey(ie,1)   0          0          0  
        0         0        0   1   ex(ie,1)   ey(ie,1)
        1  ex(ie,2) ey(ie,2)   0          0          0  
        0         0        0   1   ex(ie,2)   ey(ie,2)
        1  ex(ie,3) ey(ie,3)   0          0          0  
        0         0        0   1   ex(ie,3)   ey(ie,3)];

    B=[0 1 0 0 0 0
       0 0 0 0 0 1
       0 0 1 0 1 0]*inv(C);

    e=B*ed(i,:)';
    if colD>3 ee([1 2 4])=e; else ee=e; end

    et=[et;ee'];
    es=[es;(D*ee)'];

    ie=ie+incie;
  end
  
else
  error('Error ! Check first argument, ptype=1 or 2 allowed')
  return
end
