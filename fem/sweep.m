  function [Y]=sweep(K,C,M,p,w,b)
% Y=sweep(K,C,M,p,w,b)
%-------------------------------------------------------------
% PURPOSE
%  Compute complex frequency response function for a 
%  system of the form
%              2
%      [K+iwC-w M]y(w)=p(w)
%
% INPUT:
%    K : stiffness matrix, dim(K)= nd x nd
%    C : damping , dim(C)= nd x nd
%    M : mass matrix, dim(M)= nd x nd
%    p : load vector, dim(p)= nd x 1
%    w : circular frequency vector, dim(w)= nw x 1
%    b : boundary condition, dim(b)= npdof's
%
% OUTPUT:
%    Y : complex response matrix 
%        dim(Y)= m x nw, m : number of dof's
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom 2004-10-07   
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [nd,nd]=size(K);
  fdof=[1:nd]';
  [nw,nc]=size(w);
%
  if nargin==6
    pdof=b(:);
    fdof(pdof)=[];
    Y=zeros(nd,nw); 
    for j=1:nw;
      wj=w(j);
      Y(fdof,j)=(K(fdof,fdof)+i*wj*C(fdof,fdof)-wj*wj*M(fdof,fdof))\p(fdof);
    end
  else  
    for j=1:nw;
      wj=w(j);
      Y(:,j)=(K+i*wj*C-wj*wj*M)\p;
    end
  end
%--------------------------end--------------------------------
