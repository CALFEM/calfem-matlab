  function [s]=spectra(a,xi,dt,f)
%[s]=spectra(a,xi,dt,f)
%-------------------------------------------------------------
% PURPOSE
%  Compute seismic response spectra for elastic design.
%
% INPUT:
%    a : acceleration time history at 
%        equal time step dim(a)= ntimes x 1
%    xi: damping ratio, dim(xi)= 1 x 1
%    dt: time step dim(dt)= 1 x 1
%    f : frequency vector, dim(f)= nfreq x 1
%
% OUTPUT:
%    s : complex response matrix 
%        dim(Y)= m x nw, m : number of dof's
%-------------------------------------------------------------

% LAST MODIFIED: H Carlsson 1993-09-21   
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
[np nr]=size(f);

for i=1:np
  r=sqrt(1-xi*xi);
  w=2*pi*f(i);
  wt=w*dt;
  e=exp(-xi*wt);
  c1=2*e*cos(r*wt);
  c2=e*e;
  c3=e/r/wt*sin(r*wt);
  c4=[-c3 2*c3 -c3];
  c5=[1 -c1 c2];
  y=filter(c4,c5,a);
  s(i)=max(abs(y+a))
end  
%--------------------------end--------------------------------
