function [Ke,fe]=flw2i8e(ex,ey,ep,D,eq)
% Ke=flw2i8e(ex,ey,ep,D)
% [Ke,fe]=flw2i8e(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness (conductivity)
%  matrix for 8 node isoparametric field element
%
% INPUT:  ex = [x1 ... x8]    element coordinates
%         ey = [y1 ... y8]
%                             
%         ep = [t ir]          thickness and 
%                              integration rule
%
%         D  = [kxx kxy;
%               kyx kyy]       constitutive matrix
%
%         eq                   heat supply per unit volume
%
% OUTPUT: Ke :  element 'stiffness' matrix (8 x 8)
%         fe :  element load vector (8 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-24
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  t=ep(1); ir=ep(2); ngp=ir*ir;
  if nargin==4; eq=0 ; end

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


  Ke1=zeros(8,8);  fe1=zeros(8,1);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    JTinv=inv(JT(indx,:));
    B=JTinv*dNr(indx,:);
    Ke1=Ke1+B'*D*B*detJ*wp(i);
    fe1=fe1+N(i,:)'*detJ*wp(i);
  end

  Ke=Ke1*t;  fe=fe1*t*eq;
%--------------------------end--------------------------------
