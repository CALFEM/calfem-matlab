 function [Q,T,beta_first,beta_res] = lan(K,M,f,m)
%--------------------------------------------------------------------------
%  SYNTAX  [Q,T,beta_first,beta_res] = lan(K,M,f,m)
%
%  PURPOSE
%     Transform a pair of matrices (K,M) with size n
%     to tridiagonal form with size m
%
%  REFERENCES
%     Written: H Carlsson 1993-02-03
%     Revised: H Carlsson 1993-09-19
%
%  INPUT:
%     K : Non-singular stiffness matrix
%     M : Positive semi-definite mass matrix  
%     f : Starting vector
%     m : number of Lanczos vectors
%
%  OUTPUT:
%     Q : n x m matrix containing the Lanczos vectors
%     T : m x m Tridiagonal matrix 
%     beta_first : normalization factor for first Lanczos vector
%     beta_res : norm of residual vector 
%
%-------------------------------------------------------------------------

%
% LAST MODIFIED: H Carlsson 1993-09-19
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------------------

[nr,nc] = size(K);
Q = zeros(nr,m);
T = zeros(m,m);
%
% 0. Factorize K
[L,U] = lu(K);
%
% 1. First Lanczos vector and residual
%
y = L\f;
r = U\y;
pbar = M*r;
beta_first = sqrt(pbar'*r);
%
Q(:,1) = r/beta_first;
p=pbar/beta_first;
y = L\p;
rbar = U\y;
%
%
alfa = p'*rbar;
r = rbar - alfa*Q(:,1);
pbar = M*r;
beta = sqrt(pbar'*r);
T(1,1) = alfa;
T(2,1) = beta;
T(1,2) = beta; 
%
% 2. Lanczos vectors 2 to m
%
for j = 2:m
%
% Full reorthogonalization 
%
   for i = 1:j-2
      zold  = Q(:,i)'*p;
      Q(:,j-1) = Q(:,j-1) - zold*Q(:,i);
      z  = Q(:,i)'*pbar;
      r = r - z*Q(:,i);
   end
   z = pbar'* Q(:,j-1); 
   r = r - z*Q(:,j-1);
   p = M*Q(:,j-1);
   qnorm = p'*Q(:,j-1);
   p = p/qnorm;
   Q(:,j-1) = Q(:,j-1)/qnorm;
   T(j-1,j) =  T(j-1,j)*qnorm;
   T(j,j-1) =  T(j,j-1)*qnorm;
   z = r'*p;
   T(j-1,j-1) = T(j-1,j-1) + z ;   
   r = r - z*Q(:,j-1);
   pbar = M*r;
   beta = sqrt(pbar'*r);
%
   Q(:,j) = r/beta;
   p = pbar/beta;
   y = L\p;
   rbar = U\y;
   rroof = rbar - beta*Q(:,j-1);
   alfa = p'*rroof;
   r = rroof - alfa*Q(:,j);
   pbar = M*r;
   beta = sqrt(pbar'*r);
   T(j,j) = alfa;
   if j ~= m  
     T(j+1,j) = beta;
     T(j,j+1) = beta;
   else
     beta_res = beta;
   end
end

%----------------------- end ---------------------------------- 
