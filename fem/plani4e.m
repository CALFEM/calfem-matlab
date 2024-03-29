function [Ke,fe]=plani4e(ex,ey,ep,D,eq)
% Ke=plani4e(ex,ey,ep,D)
% [Ke,fe]=plani4e(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a 4 node isoparametric
%  element in plane strain or plane stress.
%
% INPUT:  ex = [x1 x2 x3 x4]  element coordinates
%         ey = [y1 y2 y3 y4]
%                             
%         ep =[ptype t ir]    element property 
%                               ptype: analysis type
%                               ir: integration rule
%                               t : thickness
%
%         D                   constitutive matrix
%
%         eq = [bx; by]       bx: body force in x direction
%                             by: body force in y direction
%
% OUTPUT: Ke : element stiffness matrix (8 x 8)
%         fe : equivalent nodal forces (8 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
  ptype=ep(1); t=ep(2);  ir=ep(3);  ngp=ir*ir;
  if nargin==4   b=zeros(2,1);  else  b=eq; end

%--------- gauss points --------------------------------------
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

%--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta)/4;  N(:,2)=(1+xsi).*(1-eta)/4;
  N(:,3)=(1+xsi).*(1+eta)/4;  N(:,4)=(1-xsi).*(1+eta)/4;

  dNr(1:2:r2,1)=-(1-eta)/4;     dNr(1:2:r2,2)= (1-eta)/4;
  dNr(1:2:r2,3)= (1+eta)/4;     dNr(1:2:r2,4)=-(1+eta)/4;
  dNr(2:2:r2+1,1)=-(1-xsi)/4;   dNr(2:2:r2+1,2)=-(1+xsi)/4;
  dNr(2:2:r2+1,3)= (1+xsi)/4;   dNr(2:2:r2+1,4)= (1-xsi)/4;


  Ke=zeros(8,8);
  fe=zeros(8,1);
  JT=dNr*[ex;ey]';

%--------- plane stress --------------------------------------
if ptype==1
  
  colD=size(D,2);
  if colD>3
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end
    
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:8-1)=dNx(1,:);
      B(2,2:2:8)  =dNx(2,:);
      B(3,1:2:8-1)=dNx(2,:);
      B(3,2:2:8)  =dNx(1,:);

      N2(1,1:2:8-1)=N(i,:);
      N2(2,2:2:8)  =N(i,:);

      Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
      fe=fe+N2'*b*detJ*wp(i)*t;
  end
%--------- plane strain --------------------------------------
elseif ptype==2
  
  colD=size(D,2);
  if colD>3
    Dm=D([1 2 4],[1 2 4]);
  else
    Dm=D;
  end
    
  for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:8-1)=dNx(1,:);
      B(2,2:2:8)  =dNx(2,:);
      B(3,1:2:8-1)=dNx(2,:);
      B(3,2:2:8)  =dNx(1,:);

      N2(1,1:2:8-1)=N(i,:);
      N2(2,2:2:8)  =N(i,:);

      Ke=Ke+B'*Dm*B*detJ*wp(i)*t;
      fe=fe+N2'*b*detJ*wp(i)*t;
  end
  
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
