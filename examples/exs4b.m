% example exs4b 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane truss using loops and extraction of  
%    element coordinates from a global coordinate matrix
%----------------------------------------------------------------

% REFERENCES
%     P-E Austrell 1994-03-08 
%     K-G Olsson 1995-09-28
%     O Dahlblom 2004-08-31
%----------------------------------------------------------------

echo on 

%----- Topology matrix Edof -------------------------------------

 Edof=[1   1  2  5  6;
       2   3  4  7  8;
       3   5  6  9 10;
       4   7  8 11 12;
       5   7  8  5  6;
       6  11 12  9 10;
       7   3  4  5  6;
       8   7  8  9 10;
       9   1  2  7  8;
      10   5  6 11 12];
 
%----- Stiffness matrix K and load vector f ---------------------

 K=zeros(12); 
 f=zeros(12,1);  f(11)=0.5e6*sin(pi/6);  f(12)=-0.5e6*cos(pi/6);

%----- Element properties ---------------------------------------

 A=25.0e-4;	  E=2.1e11;   ep=[E A];	

%----- Global coordinates and topology --------------------------

 Coord=[0 2;
        0 0;
        2 2;
        2 0;
        4 2;
        4 0];
 
 Dof=[ 1  2;
       3  4;
       5  6;
       7  8;
       9 10;
      11 12];
 
%----- Element coordinates -------------------------------------- 

 [Ex,Ey]=coordxtr(Edof,Coord,Dof,2);

%----- Create element stiffness matrices Ke and assemble into K -

 for i=1:10
    Ke=bar2e(Ex(i,:),Ey(i,:),ep);
    K=assem(Edof(i,:),K,Ke);
 end
 
%----- Solve the system of equations ----------------------------

 bc= [1 0;2 0;3 0;4 0];   
 [a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------

 ed=extract(Edof,a);

 for i=1:10
    N(i,:)=bar2s(ex(i,:),ey(i,:),ep,ed(i,:));
 end
 N
 
%---------------------------- end -------------------------------
 echo off

