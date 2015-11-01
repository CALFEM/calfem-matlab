function [Y]=dyna2(w2,xi,f,g,dt)
% [Y]=dyna2(w2,xi,f,g,dt)
%-------------------------------------------------------------
% PURPOSE
%  Compute dynamic solution to a set of uncoupled 
%  second-order differential equations.
%
% INPUT:
%    w2: circular frequencies squared, dim(w2)= neq x 1
%    xi: damping ratio               , dim(xi)= neq x 1
%    f : force vector                , dim(f) = neq x 1
%    g : load function in terms of straight 
%        line segments               , dim(g)=  np x 1
%    dt: output time increment  
% OUTPUT:
%    Y : dynamic solution in  
%        dim(Y)= neq x np  neq : number of eqs.
%                           np  : number of time points
%-------------------------------------------------------------

% LAST MODIFIED: H Carlsson  1993-09-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
[neq nc]=size(w2);
[np nr]=size(g)

Y=zeros(neq,np);
w=sqrt(w2);
wd=w.*sqrt(1-xi.*xi);
w3=w.^3;

ynew=Y(:,1);
ypnew=ynew;

nstep=np-1;
for j=1:nstep

  yold=ynew;
  ypold=ypnew;

  a=g(j);
  b=(g(j+1)-g(j))/dt;
 
  a0=a./w2-2*xi.*b./w3;
  a1=b./w2;
  a2=yold-a0;
  a3=(ypold+xi.*w.*a2-a1)./wd;

  ynew=a0+a1*dt+a2.*exp(-xi.*w*dt).*cos(wd*dt)+a3.*exp(-xi.*w*dt).*sin(wd*dt);
  ypnew=a1+(wd.*a3-xi.*w.*a2).*exp(-xi.*w*dt).*cos(wd*dt)-(wd.*a2+xi.*w.*a3).*exp(-xi.*w*dt).*sin(wd*dt);
  Y(:,j+1)=ynew;

end
%--------------------------end--------------------------------
