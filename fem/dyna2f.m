function [Y]=dyna2f(w2,xi,f,p,dt)
% [Y]=dyna2f(w2,xi,f,p,dt)
%-------------------------------------------------------------
% PURPOSE
%  Compute dynamic solution to a set of uncoupled 
%  second-order differential equations in the 
%  frequency domain.
%
% INPUT:
%    w2: circular frequencies squared, dim(w2)= neq x 1
%    xi: damping ratio               , dim(xi)= neq x 1
%    f : force vector                , dim(f) = neq x 1
%    p : load function in terms of force 
%        Fourier coefficients        , dim(p)=  N x 1
%    dt: output time increment  
% OUTPUT:
%    Y : dynamic solution in frequency domain  
%        dim(Y)= neq x N  neq : number of eqs.
%                          N   : number of sampled points
%-------------------------------------------------------------

% LAST MODIFIED: H Carlsson  1993-10-08
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
[neq nc]=size(w2);
N=length(p);

Y=zeros(neq,N);
w=sqrt(w2);
wbar1=2*pi/(N*dt);
r=wbar1./w;
r2=r.*r;
ii=sqrt(-1);

for i=1:neq
  for n=0:N-1
    if n<=N/2 
      Y(i,n+1)=p(n+1)*f(i)/w2(i)/(-n*n*r2(i)+2*ii*n*r(i)*xi(i)+1);
    else 
      n2=-(N-n);
      Y(i,n+1)=p(n+1)*f(i)/w2(i)/(-n2*n2*r2(i)+2*ii*n2*r(i)*xi(i)+1);     
    end
  end
end
%--------------------------end--------------------------------
