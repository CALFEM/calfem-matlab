  function [K1,f1]=statcon(K,f,cd)
% [K1,f1]=statcon(K,f,cd)
%-------------------------------------------------------------
% PURPOSE
%  Condensation of static FE-equations according to the
%  vector cd.
%
% INPUT: K : global stiffness matrix, dim(K)= nd x nd
%        f : global load vector, dim(f)= nd x 1
%
%        cd : vector containing dof's to be eliminated
%             dim(cd)= nc x 1, nc: number of condensed dof's
%
% OUTPUT: K1 : condensed stiffness matrix,
%              dim(K1)= (nd-nc) x (nd-nc)
%         f1 : condensed load vector, dim(f1)= (nd-nc) x 1
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1993-08-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [nd,nd]=size(K);
%
  aindx=[1:nd]';
  aindx(cd)=[];
  bindx=cd;
%
  Kaa=K(aindx,aindx);
  Kab=K(aindx,bindx);
  Kbb=K(bindx,bindx);
%
  fa=f(aindx);
  fb=f(bindx);
%
  K1=Kaa-Kab*inv(Kbb)*Kab';
  f1=fa-Kab*inv(Kbb)*fb;
%--------------------------end--------------------------------

