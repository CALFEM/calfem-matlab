  function [es,et,eci]=soli8s(ex,ey,ez,ep,D,ed)
% [es,et,eci]=soli8s(ex,ey,ez,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a
%  8 node (brick) isoparametric element.
%
% INPUT:  ex = [x1 x2 x3 ... x8]
%         ey = [y1 y2 y3 ... y8]  element coordinates
%         ez = [z1 z2 z3 ... z8]
%
%         ep = [Ir]               Ir: integration rule
%
%         D                       constitutive matrix
%
%         ed = [u1 u2 ..u24;      element displacement vector
%               ...........]      one row for each element
%  
% OUTPUT: es = [ sigx sigy sigz sigxy sigyz sigxz ;  
%                  ......       ...               ] 
%         element stress matrix, one row for each element  
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa   1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  ir=ep(1);  ngp=ir*ir*ir;

%--------- gauss points --------------------------------------
  if ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
  elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;

    I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]'; I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
    gp(:,1)=[I1 I1 I1]'*g1;           gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
    I1=abs(I1);  I2=abs(I2);
    w(:,1)=[I1 I1 I1]'*w1;            w(:,1)=[I2 I2 I2]'*w2+w(:,1);
    I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]'; I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
    gp(:,2)=[I1 I1 I1]'*g1;           gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
    I1=abs(I1);  I2=abs(I2);
    w(:,2)=[I1 I1 I1]'*w1;            w(:,2)=[I2 I2 I2]'*w2+w(:,2);
    I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]'; I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
    I3=abs(I1);
    gp(:,3)=[I1 I2 I3]'*g1;           gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
    w(:,3)=[I3 I2 I3]'*w1;            w(:,3)=[I2 I3 I2]'*w2+w(:,3);
  else
    disp('Used number of integration points not implemented');
    return
  end

  wp=w(:,1).*w(:,2).*w(:,3);
  xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;

%--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
  N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
  N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
  N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

  dNr(1:3:r2,1)=-(1-eta).*(1-zet);    dNr(1:3:r2,2)= (1-eta).*(1-zet);
  dNr(1:3:r2,3)= (1+eta).*(1-zet);    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
  dNr(1:3:r2,5)=-(1-eta).*(1+zet);    dNr(1:3:r2,6)= (1-eta).*(1+zet);
  dNr(1:3:r2,7)= (1+eta).*(1+zet);    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
  dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
  dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
  dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
  dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
  dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
  dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
  dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
  dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
  dNr=dNr/8.;

%--------- three dimensional case ----------------------------
  rowed=size(ed,1);
  rowex=size(ex,1);

  if rowex==1 incie=0; else incie=1; end

  es=[]; et=[]; eci=[]; ie=1;
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:);ez(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:);ez(ie,:)]';

    for i=1:ngp
      indx=[ 3*i-2; 3*i-1; 3*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:3:24-2)=dNx(1,:);
      B(2,2:3:24-1)=dNx(2,:);
      B(3,3:3:24)  =dNx(3,:);
      B(4,1:3:24-2)=dNx(2,:);
      B(4,2:3:24-1)=dNx(1,:);
      B(5,1:3:24-2)=dNx(3,:);
      B(5,3:3:24)  =dNx(1,:);
      B(6,2:3:24-1)=dNx(3,:);
      B(6,3:3:24)  =dNx(2,:);

      ee=B*ed(ied,:)';
      et=[et; ee'];
      es=[es; (D*ee)'];
    end
    
    ie=ie+incie;
 end
%--------------------------end--------------------------------
