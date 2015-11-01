function [he]=cp2s2f(ex,ey,ep)
% [he]=cp2s2f(ex,ey,ep)
%-----------------------------------------------------------
% PURPOSE
%  Compute element coupling matrix between a 4 node 
%  isoparametric acoustic element and a 2 node beam element. 
%
%     1-------------2 Beam    
%     *-------------*
%     1-------------2 Fluid
%     |             |
%     |             |
%
%
% INPUT:  ex = [x1 x2]   element coordinates
%         ey = [y1 y2]
%                             
%         ep = [t]       thickness
%
% OUTPUT: he :  element coupling matrix (6 x 2)
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

% Form local coupling matrix according, exact integration.  
% detJ = L because integration is performed for the interval [0,1].

Cle = L*[ 0         0    ;
          7/20     3/20  ;
          L/20     L/30  ;
          0        0     ;
          3/20     7/20  ;
         -L/30    -L/20] ;

% Local to global transformation matrix G 

n=b/L;

G=[n(1) n(2) 0    0    0   0;
  -n(2) n(1) 0    0    0   0;
   0    0    1    0    0   0;
   0    0    0   n(1) n(2) 0;
   0    0    0  -n(2) n(1) 0;
   0    0    0    0    0   1];

Ce=G'*Cle;    

he=t*Ce;
%-------------------------- end -------------------------------
