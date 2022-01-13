% example exs_beambar2 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a combined beam and bar structure.
%----------------------------------------------------------------

% REFERENCES
%     Ola Dahlblom 2015-11-16
%     Ola Dahlblom 2019-12-19
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%----------------------------------------------------------------
 echo on

%----- Topology -------------------------------------------------

 Edof1=[1  1  2  3  4  5  6;
        2  4  5  6  7  8  9;
        3  7  8  9 10 11 12];
 Edof2=[4 13 14  4  5;
        5 13 14  7  8];      

%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(14);	f=zeros(14,1);

%----- Element stiffness and element load matrices  -------------

 E=200e9;
 A1=4.0e-3;     A2=1.0e-3;
 I1=5.4e-5;	
 
 ep1=[E A1 I1];	 ep4=[E A2];
  
 eq1=[0 0];      eq2=[0 -10e3];

 ex1=[0 2];      ey1=[2 2];      
 ex2=[2 4];	     ey2=[2 2];
 ex3=[4 6];      ey3=[2 2];      
 ex4=[0 2];	     ey4=[0 2];
 ex5=[0 4];	     ey5=[0 2];
 
 Ke1=beam2e(ex1,ey1,ep1);
 [Ke2,fe2]=beam2e(ex2,ey2,ep1,eq2);
 [Ke3,fe3]=beam2e(ex3,ey3,ep1,eq2);
 Ke4=bar2e(ex4,ey4,ep4);
 Ke5=bar2e(ex5,ey5,ep4);
 
 %----- Assemble Ke into K ---------------------------------------

 K=assem(Edof1(1,:),K,Ke1);
 [K,f]=assem(Edof1(2,:),K,Ke2,f,fe2);
 [K,f]=assem(Edof1(3,:),K,Ke3,f,fe3);
 K=assem(Edof2(1,:),K,Ke4);
 K=assem(Edof2(2,:),K,Ke5);

%----- Solve the system of equations and compute reactions ------

 bc=[1 0; 2 0; 3 0; 13 0; 14 0];	
 [a,r]=solveq(K,f,bc)

%----- Section forces -------------------------------------------

 Ed1=extract_ed(Edof1,a);
 Ed2=extract_ed(Edof2,a);

 es1=beam2s(ex1,ey1,ep1,Ed1(1,:),eq1,11) 
 es2=beam2s(ex2,ey2,ep1,Ed1(2,:),eq2,11) 
 es3=beam2s(ex3,ey3,ep1,Ed1(3,:),eq2,11) 
 es4=bar2s(ex4,ey4,ep4,Ed2(1,:)) 
 es5=bar2s(ex5,ey5,ep4,Ed2(2,:)) 

%------------------------ end -----------------------------------
echo off