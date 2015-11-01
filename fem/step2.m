  function [Dsnap,D,V,A]=step2(K,C,M,d0,v0,ip,f,pdisp)
%Dsnap=step2x(K,C,M,d0,v0,ip,f,pdisp)
%[Dsnap,D,V,A]=step2x(K,C,M,d0,v0,ip,f,pdisp)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for dynamic solution of second-order
%  FE equations considering boundary conditions.
%
% INPUT:
%    K : stiffness matrix, dim(K)= nd x nd
%    C : damping matrix, dim(C)= nd x nd
%    M : mass matrix, dim(M)= nd x nd
%    d0 : initial vector d(0), dim(f)= nd x 1
%    v0 : initial vector v(0), dim(f)= nd x 1
%    ip : [dt tottime alfa delta [nsnap nhist t(i) ...  dof(i) ... ]] :
%         integration and output parameters
%    f : load vector, dim(f)= n x nstep+1,  
%        If dim(f)= n x 1 it is supposed that the values 
%        are kept constant during time integration.
%    pdisp : boundary condition matrix, dim(pdisp)= nbc x nstep+2,
%            where nbc = number of boundary conditions 
%            (constant or time dependent) 
%            The first column contains the degrees-of-freedom with  prescribed values
%            and the subsequent columns contain there the time history.
%            If dim(pdisp)= nbc x 2 it is supposed that the values 
%            are kept constant during time integration.
%
% OUTPUT:
%    Dsnap : displacement snapshots at 'nsnap' timesteps, specified in ip.
%          the time are also specified in ip. dim(Dsnap)=nd x nsnap
%    D :   solution matrix containing time history displacement
%          d at the selected dof's. dim(D)=nhist x nstep+1
%    V :   solution matrix containing time history of the first time derivative  
%          of d at the selected dof's. dim(V)=nhist x nstep+1
%    A :   solution matrix containing time history of the second time derivative 
%          of d at the selected dof's. dim(A)=nhist x nstep+1
%-------------------------------------------------------------

% LAST MODIFIED: G Sandberg  1994-01-29
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
% Calfem_com='step2x(K,C,M,d0,v0,ip,f,pdisp)';  Calfem_arg=[7,8];
%   errtest;   if Calfem_fail==1 return; end;
%-------------------------------------------------------------
  [nd,nd]=size(K);
  [ndc,ndc]=size(C); if (ndc==0); C=zeros(nd,nd);  end
  dt=ip(1);  tottime=ip(2);   alfa=ip(3);   delta=ip(4);
  b1 = dt*dt*0.5*(1-2*alfa);  b2 = (1-delta)*dt;
  b3 = delta*dt;              b4 = alfa*dt*dt;
  
  nstep=1;
  [nr nc]=size(f);    if (nc>1);     nstep = nc-1; end
  if nargin==7;                      bound=0;      end
  if nargin==8   
    [nr nc]=size(pdisp); if (nc>2);  nstep = nc-2; end
    bound=1;             if (nc==0); bound=0;      end
  end
  ns=tottime/dt;   if (ns < nstep | nstep==1); nstep=ns;     end

  [nr nc]=size(f);
  tf = zeros(nd,nstep+1); 
  if (nc==1);     tf(:,:)=f(:,1)*ones(1,nstep+1);  end
  if (nc>1);      tf = f;        end
 
  [nr ncip]=size(ip);
  if (ncip >= 5);
    nsnap=ip(5);             nhist=ip(6);
    lists=ip(7:6+nsnap);     listh=ip(7+nsnap:ncip);
    if (nhist > 0);     
      [nr nc]=size(listh);
      D=zeros(nc,nstep+1);   V=zeros(nc,nstep+1);   A=zeros(nc,nstep+1); 
    end
    if (nsnap > 0);          Dsnap=zeros(nd,nsnap); end
  end  
  a0=M\(tf(:,1)-C*v0-K*d0);
  if (nhist > 0);
    D(:,1) = d0(listh);   V(:,1) = v0(listh);   A(:,1) = a0(listh);
  end 
   
  tempd=zeros(nd,1);    tempv=zeros(nd,1);    tempa=zeros(nd,1);  
  fdof=[1:nd]';  
  if (bound==1);
    [nr nc]=size(pdisp); 
    if (nc==2);
      pd=pdisp(:,2)*ones(1,nstep+1);      pv=zeros(nr,nstep+1);
    end
    if (nc>2);
      pd=pdisp(:,2:nstep+2);      pv(:,1)=(pd(:,2)-pd(:,1))/dt;
%size(pd), size(pdisp),size(pv),
      pv(:,2:nstep+1)=(pd(:,2:nstep+1)-pd(:,1:nstep))/dt;
    end
    pdof=pdisp(:,1); fdof(pdof)=[];
    Keff = M(fdof,fdof)+b3*C(fdof,fdof)+b4*K(fdof,fdof);
  end
  if (bound==0);  Keff = M+b3*C+b4*K;         end 
  [L,U]=lu(Keff); 

  dnew=d0(fdof);    vnew=v0(fdof);    anew=a0(fdof); 
  isnap=1;
  for j = 1:nstep;
    time=dt*j;
    dold=dnew;      vold=vnew;      aold=anew;
    dpred=dold+dt*vold+b1*aold;     vpred=vold+b2*aold;
    if (bound==0); reff=tf(:,j+1)-C*vpred-K*dpred;  end
    if (bound==1); 
      pdeff=C(fdof,pdof)*pv(:,j+1)+K(fdof,pdof)*pd(:,j+1);       
      reff=tf(fdof,j+1)-C(fdof,fdof)*vpred-K(fdof,fdof)*dpred-pdeff;
    end
    y=L\reff;         anew=U\y;  dnew=dpred+b4*anew;  vnew=vpred+b3*anew;
    if (nhist > 0 | nsnap > 0);
      if (bound==1);  tempd(pdof)=pd(:,j+1);  tempv(pdof)=pv(:,j+1); end
      tempd(fdof)=dnew;         tempv(fdof)=vnew;         tempa(fdof)=anew;         
      if (nhist > 0);
        D(:,j+1) = tempd(listh);  V(:,j+1) = tempv(listh);  A(:,j+1) = tempa(listh); 
      end
      if (nsnap > 0  & isnap <= nsnap );
        if (time >= lists(isnap)); Dsnap(:,isnap) = tempd; isnap=isnap+1; end
      end
    end
  end
%--------------------------end--------------------------------
