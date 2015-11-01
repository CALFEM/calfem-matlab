% example exn2 
%----------------------------------------------------------------
% PURPOSE 
%    Buckling analysis of a plane frame.
%----------------------------------------------------------------

% REFERENCES
%     Susanne Heyden 95-11-01 
%     Karl-Gunnar Olsson 96-01-23
%     Ola Dahlblom 2004-09-24
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
 
 Ex=[0 0;6 6;0 6];   Ey=[4 0;4 0;4 4];

% ----- Initial loads -----

 f0=zeros(12,1);	
 f0(4)=1e3;   f0(5)=-150e3;   f0(8)=-150e3;

% ----- Increase loads until det(K)=0 -----

 j=0;
 for alpha=1:0.2:10
   j=j+1;
   N=[0.01 0 0];
   N0=[1 1 1];

% ----- Iteration for convergence -----

   eps=0.0001;		
   n=0;			
   while(abs((N(1)-N0(1))/N0(1))>eps)
     n=n+1;

     K=zeros(12,12);
     f=f0*alpha;
     Ke1=beam2g(Ex(1,:),Ey(1,:),ep1,N(1));
     Ke2=beam2g(Ex(2,:),Ey(2,:),ep1,N(2));
     Ke3=beam2g(Ex(3,:),Ey(3,:),ep3,N(3));

     K=assem(Edof(1,:),K,Ke1);
     K=assem(Edof(2,:),K,Ke2);
     K=assem(Edof(3,:),K,Ke3);

     bc=[1 0;2 0;3 0;10 0;11 0];	
     [a,r]=solveq(K,f,bc);

     Ed=extract(Edof,a);

     es1=beam2gs(Ex(1,:),Ey(1,:),ep1,Ed(1,:),N(1));
     es2=beam2gs(Ex(2,:),Ey(2,:),ep1,Ed(2,:),N(2));
     es3=beam2gs(Ex(3,:),Ey(3,:),ep3,Ed(3,:),N(3));

     N0=N;
     N=[es1(1,1),es2(1,1),es3(1,1)];
   
     if(n>20)
       disp(['Alpha= ',num2str(alpha),': The solution doesn''t converge.'])
       break
     end
   end

% ----- Check determinant for buckling -----
 
   Kred=red(K,bc(:,1));
   if (det(Kred)<=0)
     disp(['Alpha= ',num2str(alpha),': Determinant <= 0, buckling load passed.'])
     break
   end
   if(n>20)
     break
   end
   disp(['Alpha= ',num2str(alpha),' is OK! ', int2str(n), ...
   ' iterations are performed.'])
   %disp(['  '])

% ----- Save values for plotting of results -----

   deform(j)=a(4);
   M(j)=r(3);
   loadfact(j)=alpha;
   bmode=a;
 end

% ----- Plot results -----

 figure(1), clf, plot(deform(:),loadfact(:),'+',deform(:),loadfact(:),'--')
 axis([0 0.4 0 7]), grid
 xlabel('Horizontal displacement (m)'), ylabel('alpha')
 title('Displacement(alpha) for the upper left corner')

 figure(2), clf, plot(M(:),loadfact(:),'+',M(:),loadfact(:),'--')
 axis([0 0.4e6 0 7]), grid
 xlabel('Moment in A (Nm)'), ylabel('alpha')
 title('Supporting moment M-A(alpha)')

 figure(3), clf, axis off
 eldraw2(Ex,Ey,[2,3,0]);
 Ed1=extract(Edof,bmode);
 sfac=scalfact2(Ex,Ey,Ed1,0.1);
 eldisp2(Ex,Ey,Ed1,[1,1,1],sfac);
 title('Shape of buckling mode')

 echo off


