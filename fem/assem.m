function [K_new,f_new] = assem(edof,K,Ke,f,fe)
% K=assem(edof,K,Ke)
% [K,f]=assem(edof,K,Ke,f,fe)
%-------------------------------------------------------------
% PURPOSE
%  Assemble element matrices Ke ( and fe ) into the global
%  stiffness matrix K ( and the global force vector f )
%  according to the topology matrix edof.
%
% INPUT: 
%   edof:      dof topology matrix
%     Ke:      element stiffness matrix
%     K :      the global stiffness matrix
%     f :      the global force vector
%     fe:      element force vector
%
% OUTPUT: 
%   K :        the new global stiffness matrix
%   f :        the new global force vector
%-------------------------------------------------------------

% LAST MODIFIED: H Yingvar  2025-12-26 
% Change Log:
% 1. Remove the element numbers from the first column of the input parameter `edof`.
% 2. `fe` now supports matrix input.
% 3. Replacing loop operations with matrix operations significantly improves the speed of assembling large stiffness matrices.
%-------------------------------------------------------------
%% Identify input information
num_elems = size(edof,1); % number of elements to assemble
[num_K_rows,num_K_cols] = size(K); % size of the global stiffness matrix
if num_K_rows ~= num_K_cols
    error('The global stiffness matrix K must be a square matrix.');
end
[num_Ke_rows,num_Ke_cols] = size(Ke); % size of the element stiffness matrix
if num_Ke_rows ~= num_Ke_cols
    error('The element stiffness matrix Ke must be a square matrix.');
else
    num_Ke_eles = num_Ke_rows * num_Ke_cols; % number of entries in the element stiffness matrix
end

%% Assemble global stiffness matrix
iK = reshape(kron(edof,ones(num_Ke_rows,1))',num_Ke_eles*num_elems,1);
jK = reshape(kron(edof,ones(1,num_Ke_cols))',num_Ke_eles*num_elems,1);
sK = reshape(Ke(:) * ones(1,num_elems),num_Ke_eles*num_elems,1);
K_new = K + sparse(iK,jK,sK,num_K_rows,num_K_cols);

%% Assemble global force vector
if nargin == 5
    f_new = f + sparse(edof(:),1, fe(:),num_K_rows,1);
end
end