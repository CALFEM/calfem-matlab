function [coord,edof,dof]=getmesh(Geometry,eltype)
%--------------------------------------------------------
% function [coord,edof,dof]=getmesh(Geometry,eltype)
% PURPOSE
%    Creates element mesh from geometry definitions
%
% INPUT    
%    eltype     -   string containing element type 
%                      ex: 'flw2q', 'planq', 'plantc',
%                          'planr', 'plani4', 'plant'
%    Geometry   -   Geometry definitions from CalfemPre
%
% OUTPUT
%    coord      -   global coordinate matrix 
%    edof       -   topology matrix
%    dof        -   global nodal dof matrix 
%--------------------------------------------------------

% REFERENCES :   Anders K Olsson 2002-10-10
%
% Copyright (c)  Division of Structural Mechanics
%                Lund Institute of Technology
%---------------------------------------------------------

ndomains=size(Geometry.surfaces,2);
Gdof=cell(ndomains,1);

% ----- Meshing/creating domains -------------------------

Domains=gnew('domain');
for i=1:ndomains
   [Domains]=mesh2d(Geometry,Domains,i,eltype);
end;

% ----- Show element mesh --------------------------------

figure;
ddraw2(Domains);
zoom on;

% ----- Create Gdof --------------------------------------

[Gdof,Domains,dofcount]=gtopo(Geometry,Domains);

% ----- Create global edof, dof, coord -------------------

k=0; nel=0;
for i=1:ndomains
   nddof=length(Domains(i).Dof);
   [ndel,nenp]=size(Domains(i).Edof);
   dof(k+[1:nddof],:)=Domains(i).Dof;
   coord(k+[1:nddof],:)=Domains(i).Coord;
   edof(nel+[1:ndel],:)=Domains(i).Edof;
   k=k+nddof; nel=nel+ndel;
end;
edof(:,1)=[1:nel]';
nen=nenp-1;
[dof,id]=sort(dof);
coord=coord(id,:);
id=find([dof(:,1)]~=[0; dof(1:length(dof)-1,1)]);
dof=dof(id,:); coord=coord(id,:);
