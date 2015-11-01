% example exn1 
%----------------------------------------------------------------
% PURPOSE 
%    Analysis of a plane frame using second order theory.
%----------------------------------------------------------------

% REFERENCES
%     Susanne Heyden 95-11-01 
%     Karl-Gunnar Olsson 96-01-23
%     Ola Dahlblom 2004-09-23
%----------------------------------------------------------------
 echo off

% ----- Topology -----

 Edof=[1  4  5  6  1  2  3 ;
       2  7  8  9 10 11 12;
       3  4  5  6  7  8  9]; 

% ----- Element properties and global coordinates ----- 
      
 E=200e9;  
 A1=2e-3;    A2=6e-3;
 I1=1.6e-5;	 I2=5.4e-5;
 ep1=[E A1 I1];	 ep3=[E A2 I2];
 eq3=[-50e3];

 Ex=[0 0;6 6;0 6];   Ey=[4 0;4 0;4 4];

% ----- Initial values for the iteration -----

 eps=0.0001;		% Error norm
 N=[0.01 0 0];		% Initial normal forces
 N0=[1 1 1];		% Normal forces of the initial former iteration
 n=0;			% Iteration counter

% ----- Iteration procedure -----

 while(abs((N(1)-N0(1))/N0(1))>eps)
   n=n+1

   K=zeros(12,12);
   f=zeros(12,1);	
   f(4)=10e3;

   Ke1=beam2g(Ex(1,:),Ey(1,:),ep1,N(1));
   Ke2=beam2g(Ex(2,:),Ey(2,:),ep1,N(2));
   [Ke3,fe3]=beam2g(Ex(3,:),Ey(3,:),ep3,N(3),eq3);

   K=assem(Edof(1,:),K,Ke1);
   K=assem(Edof(2,:),K,Ke2);
   [K,f]=assem(Edof(3,:),K,Ke3,f,fe3);

   bc=[1 0;2 0;3 0;10 0;11 0];	
   [a,r]=solveq(K,f,bc)

   Ed=extract(Edof,a);

   es1=beam2gs(Ex(1,:),Ey(1,:),ep1,Ed(1,:),N(1))
   es2=beam2gs(Ex(2,:),Ey(2,:),ep1,Ed(2,:),N(2))
   es3=beam2gs(Ex(3,:),Ey(3,:),ep3,Ed(3,:),N(3),eq3)

   N0=N;
   N=[es1(1,1),es2(1,1),es3(1,1)];
   if(n>20)
     disp('The solution doesn''t converge')
     return
   end
 end
 echo off


