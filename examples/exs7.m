% example exs7
%----------------------------------------------------------------
% PURPOSE 
%    Set up a frame, consisting of both beams and bars, and
%    illustrate the calculations by use of graphics functions.
%----------------------------------------------------------------

% REFERENCES
%     P-A Hansson  1994-01-20
%     K-G Olsson   1995-09-28
%     O Dahlblom   2004-10-07
%----------------------------------------------------------------
 echo on 

%----- System matrices ------------------------------------------

 K=zeros(18,18);  
 f=zeros(18,1);  f(13)=1;

 Coord=[0 0;
        1 0;
        0 1;
        1 1;
        0 2;
        1 2];
       
 Dof=[1  2  3;
      4  5  6;
      7  8  9;
     10 11 12;
     13 14 15;
     16 17 18];
    
%----- Element properties, topology and coordinates -------------

 ep1=[1 1 1]; 
 Edof1=[1   1   2   3   7   8   9;
        2   7   8   9  13  14  15;
        3   4   5   6  10  11  12;
        4  10  11  12  16  17  18;
        5   7   8   9  10  11  12;
        6  13  14  15  16  17  18];
                      
 [Ex1,Ey1]=coordxtr(Edof1,Coord,Dof,2); 
        
 ep2=[1 1];  
 Edof2=[7   1   2  10  11;
        8   7   8  16  17;
        9   7   8   4   5;
       10  13  14  10  11];
                   
 [Ex2,Ey2]=coordxtr(Edof2,Coord,Dof,2);
                  
%----- Draw the fe-mesh as a check of the model -----------------
 
 eldraw2(Ex1,Ey1,[1 3 1]);     
 eldraw2(Ex2,Ey2,[1 2 1]);
  
%----- Create and assemble element matrices ---------------------

 for i=1:6
   Ke=beam2e(Ex1(i,:),Ey1(i,:),ep1);       
   K=assem(Edof1(i,:),K,Ke);
 end

 for i=1:4
    Ke=bar2e(Ex2(i,:),Ey2(i,:),ep2);        
    K=assem(Edof2(i,:),K,Ke);
 end

%----- Solve equation system ------------------------------------

 bc= [1 0; 2 0; 3 0; 4 0; 5 0; 6 0]; 
 [a,r]=solveq(K,f,bc);

%-- Extract element displacements and display the deformed mesh -

 Ed1=extract(Edof1,a);  
 Ed2=extract(Edof2,a);
 
 [sfac]=scalfact2(Ex1,Ey1,Ed1,0.1);
 eldisp2(Ex1,Ey1,Ed1,[2 1 1],sfac);
 eldisp2(Ex2,Ey2,Ed2,[2 1 1],sfac);

%-------------------------- end ---------------------------------
