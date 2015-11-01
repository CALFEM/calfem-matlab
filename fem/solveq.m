  function [d,Q]=solveq(K,f,bc)
% a=solveq(K,f)
% [a,Q]=solveq(K,f,bc)
%-------------------------------------------------------------
% PURPOSE
%  Solve static FE-equations considering boundary conditions.
%
% INPUT: K : global stiffness matrix, dim(K)= nd x nd
%        f : global load vector, dim(f)= nd x 1
%
%        bc : boundary condition matrix
%            dim(bc)= nbc x 2, nbc : number of b.c.'s
%
% OUTPUT:  a : solution including boundary values
%          Q : reaction force vector
%              dim(a)=dim(Q)= nd x 1, nd : number of dof's
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa   1993-10-06
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  if nargin==2 ; 
     d=K\f ; 
  elseif nargin==3;
     [nd,nd]=size(K);
     fdof=[1:nd]';
%
     d=zeros(size(fdof));
     Q=zeros(size(fdof));
%
     pdof=bc(:,1);
     dp=bc(:,2);
     fdof(pdof)=[];
%
     s=K(fdof,fdof)\(f(fdof)-K(fdof,pdof)*dp);
%
     d(pdof)=dp;
     d(fdof)=s;
  end  
     Q=K*d-f;
   
%--------------------------end--------------------------------
