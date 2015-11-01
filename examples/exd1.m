% example exd1
%----------------------------------------------------------------
% PURPOSE 
%    Set up the fe-model and perform eigenvalue analysis
%    for a simple frame structure.
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------
 echo on

% ------ Generate the model ------------------------------------------

% ------ material data ------------------------------------------
E=3e10;                 rho=2500;
Av=0.1030e-2;           Iv=0.0171e-4;             % IPE100
Ah=0.0764e-2;           Ih=0.00801e-4;            % IPE80
epv=[E Av Iv rho*Av];   eph=[E Ah Ih rho*Ah];
% ------ topology -----------------------------------------------
Edof=[1   1  2  3  4  5  6
      2   4  5  6  7  8  9
      3   7  8  9 10 11 12
      4  10 11 12 13 14 15];
% ------ list of coordinates  -----------------------------------
Coord=[0  0; 0  1.5; 0  3; 1  3; 2  3];
% ------ list of degrees-of-fredom  -----------------------------
Dof=[1  2  3; 4  5  6; 7  8  9; 10 11 12; 13 14 15];
% ------ generate element matrices, assemble in global matrices - 
K=zeros(15);     M=zeros(15);
[Ex,Ey]=coordxtr(Edof,Coord,Dof,2);
for i=1:2
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),epv);
  K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
end
for i=3:4
  [k,m,c]=beam2d(Ex(i,:),Ey(i,:),eph);
  K=assem(Edof(i,:),K,k);  M=assem(Edof(i,:),M,m);  
end
 

% ----- Draw a plot of the element mesh --------------------------


clf;     eldraw2(Ex,Ey,[1 2 2],Edof);
grid;    title('2-D Frame Structure') 

% ----- Eigenvalue analysis --------------------------------------

b=[1 2 3 14]';
[La,Egv]=eigen(K,M,b);
Freq=sqrt(La)/(2*pi);


% ----- plot one eigenmode ---------------------------------------

figure(1),    clf,     grid,     title('The first eigenmode'), 
eldraw2(Ex,Ey,[2 3 1]); 
Edb=extract(Edof,Egv(:,1));      eldisp2(Ex,Ey,Edb,[1 2 2]);
FreqText=num2str(Freq(1));       text(.5,1.75,FreqText);

% ----- plot eight eigenmodes ------------------------------------

figure(2), clf, axis('equal'), hold on, axis off
sfac=0.5;
title('The first eight eigenmodes (Hz)' )
for i=1:4;
  Edb=extract(Edof,Egv(:,i));
  Ext=Ex+(i-1)*3;                eldraw2(Ext,Ey,[2 3 1]); 
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  FreqText=num2str(Freq(i));     text(3*(i-1)+.5,1.5,FreqText);
end;
Eyt=Ey-4; 
for i=5:8;
  Edb=extract(Edof,Egv(:,i));
  Ext=Ex+(i-5)*3;                eldraw2(Ext,Eyt,[2 3 1]); 
  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  FreqText=num2str(Freq(i));     text(3*(i-5)+.5,-2.5,FreqText);
end

% -------------------- end ---------------------------------------
echo off
