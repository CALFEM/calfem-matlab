function [es,deps,st]=mises(ptype,mp,est,st);
% [es,deps,st]=mises(ptype,mp,est,st)
%-------------------------------------------------------------
% PURPOSE
%  Compute stresses and plastic strains for an elasto-plastic
%  isotropic hardening von Mises material.
%
% INPUT:  ptype             analysis type
%                             1:  plane stress
%                             2:  plane strain
%                             3:  axisymmery
%                             4:  three dimensional case
%
%         mp=[E,v,h]        material properties
%                             E:  modulus of elasticity
%                             v:  Poisson's ratio
%                             h:  modulus of plasticity
%
%         est               elastic trail stresses
%
%         st=[Yi,Sy,Epeff]  internal state variables
%                             Yi:    1 if yielding 0 if not
%                             Sy:    yield stress
%                             Epeff: effective plastic strain
%
% OUTPUT: es                stresses
%         deps              plastic strains
%         st=[Yi,Sy,Epeff]  updated state variables
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1995-10-15
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

G=mp(1)/2/(1+mp(2));
H=mp(3);
[nr,nc]=size(est);
ivn=st;

%---------- plane stess --------------------------------------
if ptype==1
  disp('Plane Stress not implemented');
  
%---------- plane strain and axisymmetry----------------------
elseif (ptype==2 | ptype==3)
  for i=1:nr
    Sy=ivn(i,2);
    Skk=est(i,1)+est(i,2)+est(i,3);
    S=est(i,:)-Skk*[1 1 1 0]/3;
    SeffT=sqrt(3/2*(S*S'+S(4)*S(4)));
    fT=SeffT-Sy;
    if fT>0
      Dl=(SeffT-Sy)/(3*G+H);
      Sy2=Sy+H*Dl;
      if Sy2<0   % for softening plasticity i.e. H<0
        Dl=-Sy/H+SeffT/(3*G);
        Sy2=0;
      end
      Epeff=ivn(i,3)+Dl;
      S2=S/SeffT*Sy2;
      deps(i,:)=(S-S2)/2/G;
      es(i,:)=[S2+Skk*[1 1 1 0]/3];
      st(i,:)=[1 Sy2 Epeff];
    else
      es(i,:)=est(i,:);
      deps(i,:)=zeros(1,4);
      st(i,:)=[0 ivn(i,2:3)];
    end
   end
%--------- three dimensional case ----------------------------
elseif ptype==4
 for i=1:nr
  Sy=ivn(i,2);
  Skk=est(i,1)+est(i,2)+est(i,3);
  S=est(i,:)-Skk*[1 1 1 0 0 0]/3;
  SeffT=sqrt(3/2*(S*S'+S(4)*S(4)+S(5)*S(5)+S(6)*S(6)));
  fT=SeffT-Sy;
  if fT>0
    Dl=(SeffT-Sy)/(3*G+H);
    Sy2=Sy+H*Dl;
    if Sy2<0     % for softening plasticity i.e. H<0
      Dl=-Sy/H+SeffT/(3*G);
      Sy2=0;
    end
    Epeff=ivn(i,3)+Dl;
    S2=S/SeffT*Sy2;
    deps(i,:)=(S-S2)/2/G;
    es(i,:)=[S2+Skk*[1 1 1 0 0 0]/3];
    st(i,:)=[1 Sy2 Epeff];
  else
    es(i,:)=est(i,:);
    deps(i,:)=zeros(1,4);
    st(i,:)=[0 ivn(i,2:3)];
  end
 end
else
   error('Error ! Check first argument, ptype=1,2,3 or 4 allowed')
   return
end
%--------------------------end--------------------------------

