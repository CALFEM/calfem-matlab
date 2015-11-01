function [Ke,Me,fe]=aco2i4d(ex,ey,ep,eq)
% [Ke,Me]=aco2i4d(ex,ey,ep)
% [Ke,Me,fe]=aco2i4d(ex,ey,ep,eq)
%-------------------------------------------------------------
% PURPOSE
%  Compute element stiffness and mass  
%  matrix for 4 node isoparametric acoustic element
%
% INPUT:  ex = [x1 x2 x3 x4]   element coordinates
%         ey = [y1 y2 y3 y4]
%                             
%         ep = [t c raa ir]    thickness, speed of sound, 
%                              density, integration rule
%
%         eq                   mass inflow per unit volume
%                              and time (second derivative)
%
% OUTPUT: Ke :  element 'stiffness' matrix (4 x 4)
%         Me :  element mass matrix (4 x 4)
%         fe :  element load vector (4 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: G Sandberg    1996-03-08
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  t=ep(1); c=ep(2); raa=ep(3); ir=ep(4); ngp=ir*ir;
  if nargin==3; eq=0 ; end

  if ir==1
    g1=0.0; w1=2.0;
    gp=[ g1 g1 ];  w=[ w1 w1 ];
  elseif ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-g1; g1;-g1; g1];  gp(:,2)=[-g1;-g1; g1; g1];
    w(:,1)=[ w1; w1; w1; w1];   w(:,2)=[ w1; w1; w1; w1];
  elseif ir==3
    g1=0.774596699241483; g2=0.;
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

  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;


  Ke1=zeros(4,4); Me1=zeros(4,4);  fe1=zeros(4,1);
  JT=dNr*[ex;ey]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    JTinv=inv(JT(indx,:));
    B=JTinv*dNr(indx,:);
    Ke1=Ke1+B'*B*detJ*wp(i);
    Me1=Me1+N(i,:)'*N(i,:)*detJ*wp(i);
    fe1=fe1+N(i,:)'*detJ*wp(i);
  end

  Ke=Ke1*t*c*c; Me=Me1*t;  fe=fe1*t*c*c*raa*eq;
%--------------------------end--------------------------------
