% example exs_beam1
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a simply supported beam.
%----------------------------------------------------------------

% REFERENCES
%     Ola Dahlblom 2015-11-13
%     Ola Dahlblom  2019-12-11
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%----------------------------------------------------------------
 echo on

%----- Topology -------------------------------------------------

 Edof=[1  1  2  3  4;
       2  3  4  5  6];

%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(6);   f=zeros(6,1);   f(3)=-10000;

%----- Element stiffness matrices  ------------------------------

 E=210e9;      I=2510e-8;    ep=[E I];
 ex1=[0 3];     ex2=[3 9];

 Ke1=beam1e(ex1,ep)
 Ke2=beam1e(ex2,ep)
 
%----- Assemble Ke into K ---------------------------------------

 K=assem(Edof(1,:),K,Ke1);
 K=assem(Edof(2,:),K,Ke2);

%----- Solve the system of equations and compute support forces -

 bc=[1 0; 5 0];
 [a,r]=solveq(K,f,bc)

%----- Section forces -------------------------------------------

 Ed=extract_ed(Edof,a);

 [es1,edi1]=beam1s(ex1,ep,Ed(1,:),[0],4)
 [es2,edi2]=beam1s(ex2,ep,Ed(2,:),[0],7)

 %----- Draw deformed beam ---------------------------------------
 
 figure(1)
 hold on;
 plot([0 9],[0 0]);
 c=plot([0,0:1:3,3:1:9,9],[0;edi1(:,1);edi2(:,1);0]);
 set(c,'LineWidth',[2]);
 axis([-1 10 -0.03 0.01]);
 title('displacements')
 
 %----- Draw shear force diagram----------------------------------
 
 figure(2)
 hold on;
 plot([0 9],[0 0]);
 c=plot([0,0:1:3,3:1:9,9],[0;es1(:,1);es2(:,1);0]);
 set(c,'LineWidth',[2]);
 axis([-1 10 -8000 5000]);
 set(gca, 'YDir','reverse');
 title('shear force')

 %----- Draw moment diagram----------------------------------
 
 figure(3)
 hold on;
 plot([0 9],[0 0]);
 c=plot([0,0:1:3,3:1:9,9],[0;es1(:,2);es2(:,2);0]);
 set(c,'LineWidth',[2]);
 axis([-1 10 -5000 25000]);
 set(gca, 'YDir','reverse');
 title('moment')

 %------------------------ end -----------------------------------
 echo off
