% example exd4
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, time dependent
%    boundary conditions.
%
%    Note: example exd1.m must be run first.
%
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08
%     Karl-Gunnar Olsson 1995-09-29 
%----------------------------------------------------------------
echo on

% ----- Time dependent boundary condition -----------------------

dt=0.002;         T=1;
% ----- boundary condition, initial condition -------------------
G=[0 0; 0.1 0.02; 0.2 -0.01; 0.3 0.0; T 0];   [t,g]=gfunc(G,dt);
bc=zeros(4, 1 + length(g));
bc(1,:)=[1 g];  bc(2,1)=2;   bc(3,1)=3;   bc(4,1)=14;
d0=zeros(15,1);              v0=zeros(15,1);
% ----- output parameters ---------------------------------------
ntimes=[0.1:0.1:1];          nhist=[1 4 11]; 
% ----- time integration parameters -----------------------------
ip=[dt T 0.25 0.5 10 3 ntimes nhist];
% ----- time integration ----------------------------------------
k=sparse(K);                    m=sparse(M); 
[Dsnap,D,V,A]=step2(k,[],m,d0,v0,ip,[],bc);
% ----- plot time history for two DOF:s -------------------------
figure(1), plot(t,D(1,:),'-',t,D(2,:),'--',t,D(3,:),'-.')
grid, xlabel('time (sec)'), ylabel('displacement (m)')
title('Displacement(time) at the 1st, 4th and 11th degree-of-freedom')
text(0.2,0.022,'solid line = bottom, vertical beam, x-direction')
text(0.2,0.017,'dashed line = center, vertical beam, x-direction')
text(0.2,0.012,'dashed-dotted line = center, horizontal beam, y-direction')
% ----- plot displacement for some time increments --------------
figure(2),clf, axis('equal'), hold on, axis off,  sfac=20;  
title('Snapshots (sec), magnification = 20');
for i=1:5;
  Ext=Ex+(i-1)*3;            eldraw2(Ext,Ey,[2 3 0]); 
  Edb=extract(Edof,Dsnap(:,i));
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  Time=num2str(ntimes(i));   text(3*(i-1)+.5,1.5,Time);
end;
Eyt=Ey-4; 
for i=6:10;
  Ext=Ex+(i-6)*3;            eldraw2(Ext,Eyt,[2 3 0]); 
  Edb=extract(Edof,Dsnap(:,i));
  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  Time=num2str(ntimes(i));   text(3*(i-6)+.5,-2.5,Time);
end
% ----------------------- end -----------------------------------
echo off
