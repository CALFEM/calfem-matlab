 function [ef]=plantf(ex,ey,ep,es)
% ef=plantf(ex,ey,ep,es)
%-------------------------------------------------------------
% PURPOSE
%  Compute internal element force vector in a triangular element
%  in plane stress or plane strain. 
%
% INPUT:  ex = [x1 x2 x3]         node coordinates
%         ey = [y1 y2 y3]
%
%         ep = [ptype t]          ptype: analysis type
%                                 t : thickness
%
%         es = [ sigx sigy [sigz] tauxy  element stress matrix
%                  ......              ] one row for each element
%
% OUTPUT: fe = [f1 f2 ...f8]';     internal force vector
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
ptype=ep(1);
t=ep(2);

[rowes,colD]=size(es);
rowex=size(ex,1);

%--------- plane stress --------------------------------------
if ptype==1

  if rowex==1 incie=0; else incie=1; end

  ef=[];ie=1;
  for i=1:rowes

    C=[ 1  ex(ie,1) ey(ie,1)   0          0          0  
        0         0        0   1   ex(ie,1)   ey(ie,1)
        1  ex(ie,2) ey(ie,2)   0          0          0  
        0         0        0   1   ex(ie,2)   ey(ie,2)
        1  ex(ie,3) ey(ie,3)   0          0          0  
        0         0        0   1   ex(ie,3)   ey(ie,3)];
    A=1/2*det([ones(3,1) ex(ie,:)' ey(ie,:)']);
  
    B=[0 1 0 0 0 0;
       0 0 0 0 0 1;
       0 0 1 0 1 0]*inv(C);
    
    if colD>3 stress=es(i,[1 2 4]); else stress=es(i,:); end
    ef=[ef;(A*t*B'*stress')'];
  
    ie=ie+incie;
  end
   
%--------- plane strain --------------------------------------
elseif ptype==2

  if rowex==1 incie=0; else incie=1; end

  ef=[];ie=1;
  for i=1:rowes

    C=[ 1  ex(ie,1) ey(ie,1)   0         0         0  
        0        0       0   1   ex(ie,1)   ey(ie,1)
        1  ex(ie,2) ey(ie,2)   0         0         0  
        0        0       0   1   ex(ie,2)   ey(ie,2)
        1  ex(ie,3) ey(ie,3)   0         0         0  
        0        0       0   1   ex(ie,3)   ey(ie,3)];
    A=1/2*det([ones(3,1) ex(ie,:)' ey(ie,:)']);
  
    B=[0 1 0 0 0 0
       0 0 0 0 0 1
%       0 0 0 0 0 0
       0 0 1 0 1 0]*inv(C);

    if colD>3 stress=es(i,[1 2 4]); else stress=es(i,:); end
    ef=[ef;(A*t*B'*stress')'];

    ie=ie+incie;
  end
  
else
  error('Error ! Check first argument, ptype=1 or 2 allowed')
  return
end
