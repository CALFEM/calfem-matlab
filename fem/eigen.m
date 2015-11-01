  function [L,X]=eigen(K,M,b)
% [L]=eigen(K,M)
% [L]=eigen(K,M,b)
% [L,X]=eigen(K,M)
% [L,X]=eigen(K,M,b)
%-------------------------------------------------------------
% PURPOSE
%  Solve the generalized eigenvalue problem
%  [K-LM]X = 0, considering boundary conditions.
%
% INPUT:
%    K : global stiffness matrix, dim(K)= nd x nd
%    M : global mass matrix, dim(M)= nd x nd
%    b : boundary condition matrix
%        dim(b)= nb x 1
% OUTPUT:
%    L : eigenvalues stored in a vector with length (nd-nb) 
%    X : eigenvectors dim(X)= nd x nfdof, nfdof : number of dof's
%-------------------------------------------------------------

% LAST MODIFIED: H Carlsson  1993-09-21
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [nd,nd]=size(K);
  fdof=[1:nd]';
%
  if nargin==3
    pdof=b(:);
    fdof(pdof)=[]; 
    if nargout==2
      [X1,D]=eig(K(fdof,fdof),M(fdof,fdof));
      [nfdof,nfdof]=size(X1);
      for j=1:nfdof;
        mnorm=sqrt(X1(:,j)'*M(fdof,fdof)*X1(:,j));
        X1(:,j)=X1(:,j)/mnorm;
      end
      d=diag(D);
      [L,i]=sort(d);
      X2=X1(:,i);
      X=zeros(nd,nfdof);
      X(fdof,:)=X2;
    else
      d=eig(K(fdof,fdof),M(fdof,fdof));
      L=sort(d);
    end
  else
    if nargout==2
      [X1,D]=eig(K,M);
      for j=1:nd;
        mnorm=sqrt(X1(:,j)'*M*X1(:,j));
        X1(:,j)=X1(:,j)/mnorm;
      end
      d=diag(D);
      [L,i]=sort(d);
      X=X1(:,i);
    else
      d=eig(K,M);
      L=sort(d);
    end
  end
%--------------------------end--------------------------------
