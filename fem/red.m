  function [B]=red(A,b)
% B=red(A,b)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for reducing the size of a square
%  matrix A by omitting rows and columns defined
%  by the matrix b.
%
% INPUT:
%    b : boundary condition matrix
%            dim(b)= nbc x 1, nbc : number of b's
%    A :  unreduced matrix, dim(A)= nd x nd
%
% OUTPUT:
%    B: reduced matrix
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1993-10-06
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [nd,nd]=size(A);
  fdof=[1:nd]';
%
  pdof=b(:,1);
  fdof(pdof)=[];
%
  B=A(fdof,fdof);
%--------------------------end--------------------------------
