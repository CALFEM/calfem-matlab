% example exd2
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, full system.
%
%    Note: example exd1 must be run first.
%
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%----------------------------------------------------------------
echo on

% ----- Impact, center point, vertical beam ---------------------

dt=0.002;    T=1;
% ------ the load -----------------------------------------------
G=[0 0; 0.15 1; 0.25 0; T 0];   [t,g]=gfunc(G,dt);
f=zeros(15, length(g));         f(4,:)=1000*g;
% ------ boundary condition, initial condition ------------------
bc=[1 0; 2 0; 3 0; 14 0];
d0=zeros(15,1);                 v0=zeros(15,1);
% ------ output parameters --------------------------------------
ntimes=[0.1:0.1:1];             nhist=[4 11]; 
% ------ time integration parameters ----------------------------
ip=[dt T 0.25 0.5 10 2 ntimes nhist];
% ------ time integration ---------------------------------------
k=sparse(K);                    m=sparse(M); 
[Dsnap,D,V,A]=step2(k,[],m,d0,v0,ip,f,bc);

% ----- Plot time history for two DOF:s -------------------------

figure(1), plot(t,D(1,:),'-',t,D(2,:),'--')
grid, xlabel('time (sec)'), ylabel('displacement (m)')
title('Displacement(time) at the 4th and 11th degree-of-freedom')
text(0.3,0.017,'solid line = impact point, x-direction')
text(0.3,0.012,'dashed line = center, horizontal beam, y-direction')

% ----- Plot displacements for some time increments -------------

figure(2),clf, axis('equal'), hold on, axis off
sfac=25;  
title('Snapshots (sec), magnification = 25');
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
