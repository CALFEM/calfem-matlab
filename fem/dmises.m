function D=dmises(ptype,mp,es,st);
% D=dmises(ptype,mp,es,st)
%-------------------------------------------------------------
% PURPOSE
%  Form the elasto-plastic continuum matrix for an isotropic
%  hardening von Mises material.
%
% INPUT:  ptype             analysis type
%                             1:  plane stress
%                             2:  plane strain
%                             3:  axisymmetry
%                             4:  three dimensional case
%
%         mp=[E,v,h]        material properties
%                             E:  modulus of elasticity
%                             v:  Poisson's ratio
%                             h:  modulus of plasticity
%
%         es                stresses
%
%         st=[Yi,Sy,Epeff]  internal state variables
%                             Yi:    1 if yielding 0 if not
%                             Sy:    yield stress
%                             Epeff: effective plastic strain
%
% OUTPUT: D                 D-matrix
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
iv=st;
ny=mp(2);
G=mp(1)/2/(1+ny);
H=mp(3);
[esr,esc]=size(es);
D=[];

if ptype==1
  De=hooke(ptype,mp(1),mp(2));
  if esc~=3
     error('DMISES requires es=[Sxx Syy Sxy]');
     return
  end
  
  for i=1:esr
    if iv(i,1)>0.5
      sv=es(i,:);
      Sy=iv(i,2);
      Skk=sv(1)+sv(2);
      s=sv-Skk*[1 1 0]/3;
      ss(1)=s(1)+ny*s(2);
      ss(2)=s(2)+ny*s(1)

      A=H+9*G/(2*Sy^2)*(s*s'+s(3)*s(3)+ny/(1-ny)*(sum(s.*s)+s(3)*s(3)));
      if Sy==0
         c=0;
      else
         c=9*G*G/(A*Sy^2*(1-ny^2));
      end

      Dp=c*[      ss(1)*ss(1)       ss(1)*ss(2)    (1-ny)*ss(1)*ss(3)
                  ss(2)*ss(1)       ss(2)*ss(2)    (1-ny)*ss(2)*ss(1)
            (1-ny)*s(3)*ss(1) (1-ny)*s(3)*ss(2)  (1-ny)^2*ss(3)*ss(3)];
      D=[D;De-Dp];
    else
      D=[D;De];
    end
  end
  
elseif (ptype==2  | ptype==3)
  if esc~=4
     error('DMISES requires es=[Sxx Syy Szz Sxy]');
     return
  end
     
  De=hooke(ptype,mp(1),mp(2));
  for i=1:esr
    if iv(i,1)>0.5
      sv=es(i,:);
      Sy=iv(i,2);
      Skk=sv(1)+sv(2)+sv(3);
      s=sv-Skk*[1 1 1 0]/3;

      A=H+3*G;
      if Sy==0
         c=0;
      else
         c=9*G*G/(A*Sy*Sy);
      end

      Dp=c*s'*s;
      D=[D;De-Dp];
    else
      D=[D;De];
    end
  end
  
elseif ptype==4
  De=hooke(ptype,mp(1),mp(2));
  for i=1:esr
    if iv(i,1)>0.5
      sv=es(i,:);
      Sy=iv(i,2);
      Skk=sv(1)+sv(2)+sv(3);
      s=sv-Skk*[1 1 1 0 0 0]/3;

      A=H+3*G;
      if Sy==0
         c=0;
      else
         c=9*G*G/(A*Sy*Sy);
      end

      Dp=c*s'*s;
      D=[D;De-Dp];
    else
      D=[D;De];
    end
  end
else
   error('Error ! Check first argument, ptype=1,2,3 or 4 allowed')
   return
end
%--------------------------end--------------------------------

