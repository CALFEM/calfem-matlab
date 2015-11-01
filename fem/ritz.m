  function [L,X]=ritz(K,M,f,m,b)
% [L]=ritz(K,M,f,m)
% [L]=ritz(K,M,f,m,b)
% [L,X]=ritz(K,M,f,m)
% [L,X]=ritz(K,M,f,m,b)
%-------------------------------------------------------------
% PURPOSE
%  Compute approximative eigenvalues and eigenvectors by
%  the Lanczos method.
%
% INPUT:
%    K : global stiffness matrix, dim(K)= nd x nd
%    M : global mass matrix, dim(M)= nd x nd
%    f : starting vector
%    m : number of approximative eigenvalues
%    b : boundary condition matrix dim(bc)= nb x 1
%            
% OUTPUT:
%    L : approximative eigenvalues stored in a
%        vector with length (nd-nb) 
%    X : approximative eigenvectors 
%         dim(X)= nd x nfdof, nfdof : number of free dof's
%-------------------------------------------------------------

% LAST MODIFIED: H Carlsson  1993-09-21
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%[Q,T,beta_first,beta_res] = lan(K,M,f,m)
%
  [nd,nd]=size(K);
  fdof=[1:nd]';
%
  if nargin==5
    pdof=b(:);
    fdof(pdof)=[];
    [Qm,Tm,b_1,b_res] = lan(K(fdof,fdof),M(fdof,fdof),f(fdof),m);     
    if nargout==2
      [Sm,D]=eig(Tm);
      d=diag(D);
      d=1./d;
      [L,i]=sort(d);
      X1=Qm*Sm; 
      X2=X1(:,i);
      X=zeros(nd,m);
      X(fdof,:)=X2;
    else
      d=eig(Tm);
      d=1./d;
      L=sort(d);
    end
  else
    [Qm,Tm,b_1,b_res] = lan(K,M,f,m);     
    if nargout==2
      [Sm,D]=eig(Tm);
      d=diag(D);
      d=1./d;
      [L,i]=sort(d);
      X1=Qm*Sm; 
      X=X1(:,i);
    else
      d=eig(Tm);
      d=1./d;  
      L=sort(d);
    end
  end
%--------------------------end--------------------------------
