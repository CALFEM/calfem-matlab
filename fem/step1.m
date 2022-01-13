  function [a,da,ahist,dahist]=step1(K,C,f,a0,bc,ip,times,dofs)
% [a,da]=step1(K,C,f,a0,bc,ip)
% [a,da]=step1(K,C,f,a0,bc,ip,times)
% [a,da,ahist,dahist]=step1(K,C,f,a0,bc,ip,times,dofs)
%-------------------------------------------------------------
% PURPOSE
%  Algorithm for dynamic solution of first order
%  FE equations considering boundary conditions.
%
% INPUT:
%    K : conductivity matrix, dim(K)= nd x nd
%    C : capacity matrix, dim(C)= nd x nd
%    f : load vector, dim(f)= n x nstep+1, 
%        If dim(f)= nf x 1 it is supposed that the values 
%        are kept constant during time integration.
%    a0 : initial vector a(0), dim(a0)= nd x 1
%    bc : boundary condition matrix
%         dim(bc)= nbc x nstep+2, 
%         where nbc = number of prescribed degrees-of-freedom 
%         The first column contains the degrees-of-freedom with prescribed values
%         and the subsequent columns contain there the time history.
%         If dim(bc)= nbc x 2 it is supposed that the values of the second column
%         are kept constant during time integration.
%    ip : [dt tottime alfa],  integration parameters
%    times: [t(i) ...] times at which output should be written to a and da.
%    dofs: [dof(i) ... ] dofs for which output should be written to ahist and dahist.
%
% OUTPUT:
%    a : output matrix containing values of a at all timesteps,
%        alternatively at times specified in 'times'. 
%        dim(a)=nd x nstep+1 or dim(a)=nd x ntimes
%    da : output matrix containing values of da at all timesteps,
%        alternatively at times specified in 'times'. 
%        dim(da)=nd x nstep+1 or dim(da)=nd x ntimes
%    ahist : output matrix containing time history of a at the dofs
%        specified in 'dofs'. 
%        dim(ahist)=ndofs x nstep+1
%    dahist : output matrix containing time history of da at the dofs 
%        specified in 'dofs'. 
%        dim(dahist)=ndofs x nstep+1
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2022-01-10
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  [nd,nd]=size(K);
  dt=ip(1);  
  tottime=ip(2);  
  alfa=ip(3);  
  a1=(1-alfa)*dt;  
  a2=alfa*dt;
  
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
  if (nargin<7);
     ntimes=0;
     sa=1;
     a=zeros(nd,nstep+1);
     da=zeros(nd,nstep+1);
  else
     ntimes=length(times);
     if (ntimes > 0);        
        sa=2;
        a=zeros(nd,ntimes); 
        da=zeros(nd,ntimes); 
     end
  end
  if (nargin>=8);
     ndofs=length(dofs);
     if (ndofs > 0);     
        ahist=zeros(ndofs,nstep+1); 
        dahist=zeros(ndofs,nstep+1);
     end
  else
     ndofs=0;
  end
 
  itime=1;
  
  % Calculate initial time derivative da0
  da0=C\(tf(1)-K*a0);
  % Save initial values
  if(sa==1);
     a(:,1) = a0; 
     da(:,1) = da0; 
  elseif(sa==2);
     if (times(itime)== 0);  
        a(:,itime) = a0; 
        da(:,itime) = da0; 
        itime=itime+1; 
     end
  end 
  if (ndofs > 0);  
     ahist(:,1) = a0(dofs); 
     dahist(:,1) = da0(dofs);  
  end 
  % Reduce matrices due to bcs  
  tempa=zeros(nd,1);   
  tempda=zeros(nd,1);
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
     Keff=C(fdof,fdof)+a2*K(fdof,fdof);
  end
  if (bound==0);    
      Keff=C+a2*K;       
  end 
  [L,U]=lu(Keff);
  anew=a0(fdof);  
  danew=da0(fdof);

  % Iterate over time steps
  for j = 1:nstep
     time=dt*j;
     aold=anew;    
     daold=danew;    
     apred=aold+a1*daold;
     if (bound==0);    
        reff=tf(:,j+1)-K*apred;  
     end
     if (bound==1); 
        pdeff=C(fdof,pdof)*pda(:,j+1)+K(fdof,pdof)*pa(:,j+1);
        reff = tf(fdof,j+1)-K(fdof,fdof)*apred-pdeff;
     end
     y=L\reff;         
     danew=U\y;      
     anew=apred+a2*danew;
    % Save to a, da, ahist and dahist
     if (bound==1);  
        tempa(pdof)=pa(:,j+1);  
        tempda(pdof)=pda(:,j+1); 
     end
     tempa(fdof)=anew;         
     tempda(fdof)=danew;
     if(sa==1);
        a(:,j+1) = tempa;  
        da(:,j+1) = tempda; 
     elseif(sa==2);
        if ((ntimes > 0)  && (itime <= ntimes));
           if (time >= times(itime)); 
              a(:,itime) = tempa;  
              da(:,itime) = tempda; 
              itime=itime+1; 
           end
        end
        if (ndofs > 0);
           ahist(:,j+1) = tempa(dofs);  
           dahist(:,j+1) = tempda(dofs);
        end
     end
  end
%--------------------------end--------------------------------
