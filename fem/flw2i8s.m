function [es,et,eci]=flw2i8s(ex,ey,ep,D,ed)
% [es,et,eci]=flw2i8s(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Compute flows or corresponding quantities in the
%  8 node isoparametric field element.
%
% INPUT:  ex = [x1 x2 x3.... x8]   element coordinates
%         ey = [y1 y2 y3.... y8]
%                             
%         ep =[t ir]               thickness and 
%                                  integration rule
%
%         D = [kxx kxy;
%              kyx kyy]            constitutive matrix
%
%         ed = [u1 .... u8]        u1...u8: nodal values
%
% OUTPUT: es=[qx qy ;
%             ... ..]               element flows
%
%         et=[gx gy ;
%             ... ..]               element gradients
%
%         eci=[ix1   iy1;            Gauss point  location vector
%             ....     ;            nint: number of integration
%             ix(nint) iy(nint)]          points
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-24
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  t=ep(1); ir=ep(2); ngp=ir*ir;

  if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
  elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
  elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;
    gp(:,1)=[-g1;-g2; g1;-g1; g2; g1;-g1; g2; g1];
    gp(:,2)=[-g1;-g1;-g1; g2; g2; g2; g1; g1; g1];
    w(:,1)=[ w1; w2; w1; w1; w2; w1; w1; w2; w1];
    w(:,2)=[ w1; w1; w1; w2; w2; w2; w1; w1; w1];
  else
    disp('Used number of integration points not implemented');
    return
  end
  wp=w(:,1).*w(:,2);

  xsi=gp(:,1);  eta=gp(:,2);  r2=ngp*2;

  N(:,1)=-(1-xsi).*(1-eta).*(1+xsi+eta)/4; N(:,5)=(1-xsi.*xsi).*(1-eta)/2;
  N(:,2)=-(1+xsi).*(1-eta).*(1-xsi+eta)/4; N(:,6)=(1+xsi).*(1-eta.*eta)/2;
  N(:,3)=-(1+xsi).*(1+eta).*(1-xsi-eta)/4; N(:,7)=(1-xsi.*xsi).*(1+eta)/2;
  N(:,4)=-(1-xsi).*(1+eta).*(1+xsi-eta)/4; N(:,8)=(1-xsi).*(1-eta.*eta)/2;

  dNr(1:2:r2,1)=-(-(1-eta).*(1+xsi+eta)+(1-xsi).*(1-eta))/4;
  dNr(1:2:r2,2)=-( (1-eta).*(1-xsi+eta)-(1+xsi).*(1-eta))/4;
  dNr(1:2:r2,3)=-( (1+eta).*(1-xsi-eta)-(1+xsi).*(1+eta))/4;
  dNr(1:2:r2,4)=-(-(1+eta).*(1+xsi-eta)+(1-xsi).*(1+eta))/4;
  dNr(1:2:r2,5)=-xsi.*(1-eta);
  dNr(1:2:r2,6)=(1-eta.*eta)/2;
  dNr(1:2:r2,7)=-xsi.*(1+eta);
  dNr(1:2:r2,8)=-(1-eta.*eta)/2;
  dNr(2:2:r2+1,1)=-(-(1-xsi).*(1+xsi+eta)+(1-xsi).*(1-eta))/4;
  dNr(2:2:r2+1,2)=-(-(1+xsi).*(1-xsi+eta)+(1+xsi).*(1-eta))/4;
  dNr(2:2:r2+1,3)=-( (1+xsi).*(1-xsi-eta)-(1+xsi).*(1+eta))/4;
  dNr(2:2:r2+1,4)=-( (1-xsi).*(1+xsi-eta)-(1-xsi).*(1+eta))/4;
  dNr(2:2:r2+1,5)=-(1-xsi.*xsi)/2;
  dNr(2:2:r2+1,6)=-eta.*(1+xsi);
  dNr(2:2:r2+1,7)=(1-xsi.*xsi)/2;
  dNr(2:2:r2+1,8)=-eta.*(1-xsi);


  eci=N*[ex;ey]';  [red,ced]=size(ed);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    JTinv=inv(JT(indx,:));
    B=JTinv*dNr(indx,:);
    p1=-D*B*ed';
    p2=B*ed';
    es(i:ngp:ngp*red,:)=p1';
    et(i:ngp:ngp*red,:)=p2';
  end
%--------------------------end--------------------------------
