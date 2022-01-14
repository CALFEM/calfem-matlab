% example exd_beam2_t
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, full system.
%
%    Note: example exd_beam2_m must be run first.
%
%----------------------------------------------------------------

% REFERENCES
%     Göran Sandberg 1994-03-08 
%     Karl-Gunnar Olsson 1995-09-29
%     Ola Dahlblom 2022-01-13
%----------------------------------------------------------------
echo on

% ----- Impact, center point, vertical beam ---------------------

dt=0.002;    T=1;
% ------ the load -----------------------------------------------
G=[0 0; 0.15 1; 0.25 0; T 0];   [t,g]=gfunc(G,dt);
f=zeros(15, length(g));         f(4,:)=1000*g;
% ------ boundary condition, initial condition ------------------
bc=[1 0; 2 0; 3 0; 14 0];
a0=zeros(15,1);                 da0=zeros(15,1);
% ------ output parameters --------------------------------------
times=[0.1:0.1:1];             dofs=[4 11]; 
% ------ time integration parameters ----------------------------
ip=[dt T 0.25 0.5];
% ------ time integration ---------------------------------------
k=sparse(K);                    m=sparse(M); 
[a,da,d2a,ahist,dahist,d2ahist]=step2(k,[],m,f,a0,da0,bc,ip,times,dofs);

% ----- Plot time history for two DOF:s -------------------------

figure(1), plot(t,ahist(1,:),'-',t,ahist(2,:),'--')
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
  Edb=extract_ed(Edof,a(:,i));
  eldisp2(Ext,Ey,Edb,[1 2 2],sfac);
  Time=num2str(times(i));   text(3*(i-1)+.5,1.5,Time);
end;
Eyt=Ey-4; 
for i=6:10;
  Ext=Ex+(i-6)*3;            eldraw2(Ext,Eyt,[2 3 0]); 
  Edb=extract_ed(Edof,a(:,i));
  eldisp2(Ext,Eyt,Edb,[1 2 2],sfac);
  Time=num2str(times(i));   text(3*(i-6)+.5,-2.5,Time);
end

% ----------------------- end -----------------------------------
echo off
