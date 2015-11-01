  function [Tsnap,D,V]=step1(K,C,d0,ip,f,pbound)
% Tsnap=step1(K,C,d0,ip,f,bc)
% [Tsnap,D,V]=step1(K,C,d0,ip,f,bc)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for dynamic solution of first order
%  FE equations considering boundary conditions.
%
% INPUT:
%    K : conductivity matrix, dim(K)= nd x nd
%    C : capacity matrix, dim(C)= nd x nd
%    d0 : initial vector d(0), dim(d0)= nd x 1
%    ip : [dt tottime alfa [nsnap nhist  t(i) ...  dof(i) ... ]] :
%         integration and output parameters
%    f : load vector, dim(f)= n x nstep+1, 
%        If dim(f)= nf x 1 it is supposed that the values 
%        are kept constant during time integration.
%    pbound : boundary condition matrix
%        dim(pbound)= nbc x nstep+2, 
%        where nbc = number of prescribed degrees-of-freedom 
%        the first column contains the degrees-of-freedom with  prescribed values
%        and the subsequent columns contain there the time history.
%        If dim(pbound)= nbc x 2 it is supposed that the values of the second column
%        are kept constant during time integration.
%
% OUTPUT:
%    Tsnap : output matrix containing snapshots at 'nsnap' timesteps,
%        specified in ip, the time are also specified in ip. 
%        dim(Tsnap)=nd x nsnap
%    D : solution matrix containing time history of d at the selected dof's
%        specified in ip. dim(D)=nhist x nstep+1
%    V : solution matrix containing time history of the first time 
%        derivative d at the selected dof's specified in ip. 
%        dim(V)=nhist x nstep+1
%-------------------------------------------------------------

% LAST MODIFIED: G Sandberg  1994-01-29
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  [nd,nd]=size(K);
  dt=ip(1);  tottime=ip(2);  alfa=ip(3);  
  a1=(1-alfa)*dt;  a2=alfa*dt;
  
  nstep=1;
  [nr nc]=size(f); if (nc>1);         nstep = nc-1; end
  if nargin==5;                       bound=0;      end
  if nargin==6;   
    [nr nc]=size(pbound); if (nc>2);  nstep = nc-2; end
    bound=1;              if (nc==0); bound=0;      end
  end
  ns=tottime/dt;   if (ns < nstep | nstep==1); nstep=ns;     end
  
  [nr nc]=size(f);
  tf = zeros(nd,nstep+1); 
  if (nc==1);     tf(:,:)=f(:,1)*ones(1,nstep+1);  end
  if (nc>1);      tf = f;        end
  
  [nr ncip]=size(ip);
  if (ncip >= 4);
    nsnap=ip(4);             nhist=ip(5);
    lists=ip(6:5+nsnap);     listh=ip(6+nsnap:ncip);
    if (nhist > 0);     
      [nr nc]=size(listh); D=zeros(nc,nstep+1); V=zeros(nc,nstep+1);
    end
    if (nsnap > 0);        Tsnap=zeros(nd,nsnap); end
  end  
  
  v0=C\(tf(1)-K*d0);
  if (nhist > 0);  D(:,1) = d0(listh); V(:,1) = v0(listh);  end 
  tempd=zeros(nd,1);   tempv=zeros(nd,1);
  fdof=[1:nd]';
  if (bound==1);
    [nr nc]=size(pbound); 
    if (nc==2);
      pd=pbound(:,2)*ones(1,nstep+1);      pv=zeros(nr,nstep+1);
    end
    if (nc>2);
      pd=pbound(:,2:nstep+2);      pv(:,1)=(pd(:,2)-pd(:,1))/dt;
      pv(:,2:nstep+1)=(pd(:,2:nstep+1)-pd(:,1:nstep))/dt;
    end
    pdof=pbound(:,1); fdof(pdof)=[];
    Keff=C(fdof,fdof)+a2*K(fdof,fdof);
  end
  if (bound==0);    Keff=C+a2*K;       end 
  [L,U]=lu(Keff);
  dnew=d0(fdof);  vnew=v0(fdof);

  isnap=1;
  for j = 1:nstep
    time=dt*j;
    dold=dnew;    vold=vnew;    dpred=dold+a1*vold;
    if (bound==0);    reff=tf(:,j+1)-K*dpred;  end
    if (bound==1); 
      pdeff=C(fdof,pdof)*pv(:,j+1)+K(fdof,pdof)*pd(:,j+1);
      reff = tf(fdof,j+1)-K(fdof,fdof)*dpred-pdeff;
    end
    y=L\reff;         vnew=U\y;      dnew=dpred+a2*vnew;
    if (nhist > 0 | nsnap > 0);
      if (bound==1);  tempd(pdof)=pd(:,j+1);  tempv(pdof)=pv(:,j+1); end
      tempd(fdof)=dnew;         tempv(fdof)=vnew;
      if (nhist > 0);
        D(:,j+1) = tempd(listh);  V(:,j+1) = tempv(listh);
      end
      if ((nsnap > 0)  & (isnap <= nsnap));
        if (time >= lists(isnap)); Tsnap(:,isnap) = tempd; isnap=isnap+1; end
      end
    end
  end
end
%--------------------------end--------------------------------
