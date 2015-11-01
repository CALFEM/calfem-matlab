 function [K,f]=assem(edof,K,Ke,f,fe)
% K=assem(edof,K,Ke)
% [K,f]=assem(edof,K,Ke,f,fe)
%-------------------------------------------------------------
% PURPOSE
%  Assemble element matrices Ke ( and fe ) into the global
%  stiffness matrix K ( and the global force vector f )
%  according to the topology matrix edof.
%
% INPUT: edof:       dof topology matrix
%        K :         the global stiffness matrix
%        Ke:         element stiffness matrix
%        f :         the global force vector
%        fe:         element force vector
%
% OUTPUT: K :        the new global stiffness matrix
%         f :        the new global force vector
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa   1993-10-06
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
    [nie,n]=size(edof);
    t=edof(:,2:n);
    for i = 1:nie
      K(t(i,:),t(i,:)) = K(t(i,:),t(i,:))+Ke;
      if nargin==5
         f(t(i,:))=f(t(i,:))+fe;
      end
    end
%--------------------------end--------------------------------

