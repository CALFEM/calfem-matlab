  function [a,da,d2a,ahist,dahist,d2ahist]=step2(K,C,M,f,a0,da0,bc,ip,times,dofs)
% [a,da,d2a]=step2(K,C,M,f,a0,da0,bc,ip)
% [a,da,d2a]=step2(K,C,M,f,a0,da0,bc,ip,times)
% [a,da,d2a,ahist,dahist,d2ahist]=step2(K,C,M,f,a0,da0,bc,ip,times,dofs)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for dynamic solution of second-order
%  FE equations considering boundary conditions.
%
% INPUT:
%    K : stiffness matrix, dim(K)= nd x nd
%    C : damping matrix, dim(C)= nd x nd
%    M : mass matrix, dim(M)= nd x nd
%    f : load vector, dim(f)= n x nstep+1,  
%        If dim(f)= n x 1 it is supposed that the values 
%        are kept constant during time integration.
%    a0 : initial vector a(0), dim(a0)= nd x 1
%    da0 : initial vector v(0), dim(da0)= nd x 1
%    bc : boundary condition matrix, dim(bc)= nbc x nstep+2,
%            where nbc = number of prescribed degrees-of-freedom 
%            (constant or time dependent) 
%            The first column contains the degrees-of-freedom with prescribed values
%            and the subsequent columns contain there the time history.
%            If dim(bc)= nbc x 2 it is supposed that the values of the second column
%            are kept constant during time integration.
%    ip : [dt tottime alfa],  integration parameters
%    times: [t(i) ...] times at which output should be written to a, da and d2a.
%    dofs: [dof(i) ... ] dofs for which output should be written to ahist, dahist and d2ahist.
%
% OUTPUT:
%    a : output matrix containing values of a at all timesteps,
%        alternatively at times specified in 'times'. 
%        dim(a)=nd x nstep+1 or dim(a)=nd x ntimes
%    da : output matrix containing values of da at all timesteps,
%        alternatively at times specified in 'times'. 
%        dim(da)=nd x nstep+1 or dim(da)=nd x ntimes
%    d2a : output matrix containing values of d2a at all timesteps,
%        alternatively at times specified in 'times'. 
%        dim(d2a)=nd x nstep+1 or dim(d2a)=nd x ntimes
%    ahist : output matrix containing time history of a at the dofs
%        specified in 'dofs'. 
%        dim(ahist)=ndofs x nstep+1
%    dahist : output matrix containing time history of da at the dofs 
%        specified in 'dofs'. 
%        dim(dahist)=ndofs x nstep+1
%    d2ahist : output matrix containing time history of d2a at the dofs 
%        specified in 'dofs'. 
%        dim(dahist)=ndofs x nstep+1
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2022-01-10
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  [nd,nd]=size(K);
  [ndc,ndc]=size(C); 
  if (ndc==0); 
     C=zeros(nd,nd);  
  end
  dt=ip(1);  
  tottime=ip(2);   
  alfa=ip(3);   
  delta=ip(4);
  b1 = dt*dt*0.5*(1-2*alfa);  
  b2 = (1-delta)*dt;
  b3 = delta*dt;              
  b4 = alfa*dt*dt;
  
  nstep=1;
  [nr nc]=size(f);    
  if (nc>1);     
     nstep = nc-1; 
  end
  
  [nr nc]=size(bc); 
  if (nc>2);  
     nstep = nc-2; 
  end
  bound=1;             
  if (nc==0); 
     bound=0;      
  end
  
  ns=tottime/dt;   
  if (ns < nstep | nstep==1); 
     nstep=ns;     
  end

  [nr nc]=size(f);
  tf = zeros(nd,nstep+1); 
  if (nc==1);     
     tf(:,:)=f(:,1)*ones(1,nstep+1);  
  end
  if (nc>1);      
     tf = f;        
  end
 
  sa=0;
  if (nargin<9);
     ntimes=0;
     sa=1;
     a=zeros(nd,nstep+1);
     da=zeros(nd,nstep+1);
     d2a=zeros(nd,nstep+1);
  else
     ntimes=length(times);
     if (ntimes > 0);        
        sa=2;
        a=zeros(nd,ntimes); 
        da=zeros(nd,ntimes); 
        d2a=zeros(nd,ntimes); 
     end
  end   
  if (nargin>=10);
     ndofs=length(dofs);
     if (ndofs > 0);     
        ahist=zeros(ndofs,nstep+1); 
        dahist=zeros(ndofs,nstep+1);
        d2ahist=zeros(ndofs,nstep+1);
     end
  else
       ndofs=0;
  end
 
  itime=1;
  
  % Calculate initial second time derivative d2a0
  d2a0=M\(tf(:,1)-C*da0-K*a0);
  % Save initial values
  if(sa==1);
     a(:,1) = a0; 
     da(:,1) = da0; 
  elseif(sa==2);
     if (times(itime)== 0);  
        a(:,itime) = a0; 
        da(:,itime) = da0; 
        d2a(:,itime) = d2a0; 
        itime=itime+1; 
     end
  end   
  if (ndofs > 0);
     ahist(:,1) = a0(dofs);   
     dahist(:,1) = da0(dofs);   
     d2ahist(:,1) = d2a0(dofs);
  end 
  % Reduce matrices due to bcs  
  tempa=zeros(nd,1);    
  tempda=zeros(nd,1);    
  tempd2a=zeros(nd,1);  
  fdof=[1:nd]';  
  if (bound==1);
     [nr nc]=size(bc); 
     if (nc==2);
        pa=bc(:,2)*ones(1,nstep+1);      
        pda=zeros(nr,nstep+1);
     end
     if (nc>2);
        pa=bc(:,2:nstep+2);      
        pda(:,1)=(pa(:,2)-pa(:,1))/dt;
        pda(:,2:nstep+1)=(pa(:,2:nstep+1)-pa(:,1:nstep))/dt;
     end
     pdof=bc(:,1); fdof(pdof)=[];
     Keff = M(fdof,fdof)+b3*C(fdof,fdof)+b4*K(fdof,fdof);
  end
  if (bound==0);  
     Keff = M+b3*C+b4*K;         
  end 
  [L,U]=lu(Keff); 
  anew=a0(fdof);    
  danew=da0(fdof);    
  d2anew=d2a0(fdof); 
  
  % Iterate over time steps
  for j = 1:nstep;
     time=dt*j;
     aold=anew;      
     daold=danew;      
     d2aold=d2anew;
     apred=aold+dt*daold+b1*d2aold;     
     dapred=daold+b2*d2aold;
     if (bound==0); 
        reff=tf(:,j+1)-C*dapred-K*apred;  
     end
     if (bound==1); 
        pdeff=C(fdof,pdof)*pda(:,j+1)+K(fdof,pdof)*pa(:,j+1);       
        reff=tf(fdof,j+1)-C(fdof,fdof)*dapred-K(fdof,fdof)*apred-pdeff;
     end
     y=L\reff;         
     d2anew=U\y;  
     anew=apred+b4*d2anew;  
     danew=dapred+b3*d2anew;
    % Save to a, da, d2a, ahist, dahist and d2ahist
    if (bound==1);  
       tempa(pdof)=pa(:,j+1);  
       tempda(pdof)=pda(:,j+1); 
    end
    tempa(fdof)=anew;         
    tempda(fdof)=danew;         
    tempd2a(fdof)=d2anew;         
    if(sa==1);
       a(:,j+1) = tempa;  
       da(:,j+1) = tempda; 
       d2a(:,j+1) = tempd2a; 
    elseif(sa==2);
       if (ntimes > 0 && itime <= ntimes );
          if (time >= times(itime)); 
             a(:,itime) = tempa; 
             da(:,itime) = tempda; 
             d2a(:,itime) = tempd2a; 
             itime=itime+1;
          end   
       end
       if (ndofs > 0);
          ahist(:,j+1) = tempa(dofs);  
          dahist(:,j+1) = tempda(dofs);  
          d2ahist(:,j+1) = tempd2a(dofs); 
       end
    end
 end
%--------------------------end--------------------------------
