% example exs1 
%----------------------------------------------------------------
% PURPOSE 
%    Linear elastic spring analysis. Introduction to the basic 
%    steps in the finite element method.
%----------------------------------------------------------------

% REFERENCES
%     P-E Austrell 1994-03-08 
%     K-G Olsson 1995-09-28
%     O Dahlblom 2004-09-06
%----------------------------------------------------------------
echo on
 
%----- Topology matrix Edof -------------------------------------

 Edof=[1 1 2;
       2 2 3;
       3 2 3];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(3,3) 
 f=zeros(3,1);  f(2)=100
 
%----- Element stiffness matrices  ------------------------------

 k=1500;  ep1=k;  ep2=2*k;
 Ke1=spring1e(ep1)
 Ke2=spring1e(ep2)
 
%----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke2)
 K=assem(Edof(2,:),K,Ke1) 
 K=assem(Edof(3,:),K,Ke2)
 
%----- Solve the system of equations ----------------------------

 bc= [1 0; 3 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------

 ed1=extract(Edof(1,:),a)
 ed2=extract(Edof(2,:),a)
 ed3=extract(Edof(3,:),a)

 es1=spring1s(ep2,ed1)
 es2=spring1s(ep1,ed2)
 es3=spring1s(ep2,ed3)

%---------------------------- end -------------------------------
echo off
