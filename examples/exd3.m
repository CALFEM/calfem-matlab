% example exd3
%----------------------------------------------------------------
% PURPOSE 
%    Structural Dynamics, time integration, reduced system.
%
%    Note: example exd1.m must be run first.
%
%----------------------------------------------------------------

% REFERENCES
%     G"oran Sandberg 1994-03-08
%     Karl-Gunnar Olsson 1995-09-29 
%----------------------------------------------------------------
figure(1); clf; figure(2); clf;
echo on

% ----- Impact, center point, vertical beam ---------------------
dt=0.002;      T=1;      nev=2;
% ----- the load ------------------------------------------------
G=[0 0; 0.15 1; 0.25 0; T 0];        [t,g]=gfunc(G,dt);
f=zeros(15, length(g));              f(4,:)=1000*g;
fr=sparse([[1:1:nev]' Egv(:,1:nev)'*f]);
% ----- reduced system matrices ---------------------------------
kr=sparse(diag(diag(Egv(:,1:nev)'*K*Egv(:,1:nev))));
mr=sparse(diag(diag(Egv(:,1:nev)'*M*Egv(:,1:nev))));
% ----- initial condition ---------------------------------------
dr0=zeros(nev,1);                    vr0=zeros(nev,1);
% ----- output parameters ---------------------------------------
ntimes=[0.1:0.1:1];    nhistr=[1:1:nev];   nhist=[4 11];
% ----- time integration parameters -----------------------------
ip=[dt T 0.25 0.5 10 nev ntimes nhistr];
% ----- time integration ----------------------------------------
[Dsnapr,Dr,Vr,Ar]=step2(kr,[],mr,dr0,vr0,ip,fr,[]);
% ----- mapping back to original coordinate system --------------
DsnapR=Egv(:,1:nev)*Dsnapr;
DR=Egv(nhist,1:nev)*Dr;
% ----- plot time history for two DOF:s -------------------------
figure(1), plot(t,DR(1,:),'-',t,DR(2,:),'--')
axis([0    1.0000   -0.0100    0.0200])
grid, xlabel('time (sec)'), ylabel('displacement (m)')
title('Displacement(time) at the 4th and 11th degree-of-freedom')
text(0.3,0.017,'solid line = impact point, x-direction')
text(0.3,0.012,'dashed line = center, horizontal beam, y-direction')
text(0.3,-0.007,'TWO EIGENVECTORS ARE USED')

% ---------------------- end ------------------------------------
echo off
