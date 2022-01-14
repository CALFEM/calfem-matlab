function [K,f]=sparse_assem(edof,K,Ke,f,fe)
% K=sparse_assem(edof,K,Ke)
% [K,f]=sparse_assem(edof,K,Ke,f,fe)
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

% LAST MODIFIED: E Serrano 1997-10-10
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
ndof=size(K,1);
Ke=sparse(Ke);
[r,c]=size(edof);
rowi=zeros(r*(c-1)^2,1);
coli=zeros(r*(c-1)^2,1);

for elem=1:r
     row_p=ones(c-1,1)*edof(elem,2:c);
     rowi((elem-1)*(c-1)^2+1:elem*(c-1)^2,1)=reshape(row_p,(c-1)^2,1);
     coli((elem-1)*(c-1)^2+1:elem*(c-1)^2,1)=reshape(row_p',(c-1)^2,1);
end

%Ketmp=reshape([ones(r,1)*reshape(Ke,1,(c-1)^2)]',r*(c-1)^2,1);
%fetmp=reshape([ones(r,1)*fe']',r*(c-1),1);

Ketmp2=sparse(rowi,coli,reshape([ones(r,1)*reshape(Ke',1,(c-1)^2)]',r*(c-1)^2,1),ndof,ndof);
if nargin==5
  fetmp2=sparse(reshape(edof(:,2:c)',(c-1)*r,1),1,reshape([ones(r,1)*fe']',r*(c-1),1),ndof,1);
end;
  %[nie,n]=size(edof);
    %t=edof(:,2:n);
    %for i = 1:nie
    % Ktmp(t(i,:),t(i,:)) = Ktmp(t(i,:),t(i,:))+Ke;
    %  if nargin==5
    %      ftmp(t(i,:))=ftmp(t(i,:))+fe;
    %    end
    %   end
    %K=K+Ktmp;
K=K+Ketmp2;
 if nargin==5
    % f=f+ftmp;
  f=f+fetmp2;
 end;
%--------------------------end--------------------------------

