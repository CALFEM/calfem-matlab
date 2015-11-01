 function [es,et,eci]=plani8s(ex,ey,ep,D,ed)
% [es,et,eci]=plani8s(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear stress for a 8 node 
%  isoparametric element in plane strain or plane stress.
%
% INPUT:   ex = [x1 ... x8]    element coordinates
%          ey = [y1 ... y8]
%              
%          ep = [ptype t It ]  ptype: analysis type
%                              t : thickness 
%                              It: integration rule
%                             
%          D                   constitutive matrix
%
%          ed = [u1 u2 ..u16;   element displacement vector
%               ..........]    one row for each element
%
% OUTPUT: es = [ sigx sigy [sigz] tauxy    element stress matrix
%                  ......              ]   one row for each element
%
%         et = [ epsx epsy [epsz] gamxy    element strain matrix
%                  ......              ]   one row for each element
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  ptype=ep(1); t=ep(2);  ir=ep(3);  ngp=ir*ir;

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
  
%--------- plane stress --------------------------------------
if ptype==1
  
  rowed=size(ed,1);
  rowex=size(ex,1);
  colD =size(D ,2);
  
  if colD>3
    Cm=inv(D);
    Dm=inv(Cm([1 2 4],[1 2 4]));
  else
    Dm=D;
  end
  
  if rowex==1 incie=0; else incie=1; end
  
  es=[]; et=[]; eci=[]; ie=1;
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:16-1)=dNx(1,:);
      B(2,2:2:16)  =dNx(2,:);
      B(3,1:2:16-1)=dNx(2,:);
      B(3,2:2:16)  =dNx(1,:);

      ee=B*ed(ied,:)';
      if colD>3
         ss=zeros(colD,1);
         ss([1 2 4])=Dm*ee;
         ee=Cm*ss;
      else
         ss=Dm*ee;
      end
      
      et=[et; ee'];
      es=[es; ss'];
    end
    
    ie=ie+incie;
  end
%--------- plane stress --------------------------------------
elseif ptype==2
  
  rowed=size(ed,1);
  rowex=size(ex,1);
  colD =size(D ,2);
    
  if rowex==1 incie=0; else incie=1; end
  
  es=[]; et=[]; eci=[]; ie=1; ee=zeros(colD,1);
  for ied=1:rowed
    eci=[eci N*[ex(ie,:);ey(ie,:)]']; 
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    for i=1:ngp
      indx=[ 2*i-1; 2*i ];
      detJ=det(JT(indx,:));
      if detJ<10*eps
        disp('Jacobideterminant equal or less than zero!')
      end
      JTinv=inv(JT(indx,:));
      dNx=JTinv*dNr(indx,:);

      B(1,1:2:16-1)=dNx(1,:);
      B(2,2:2:16)  =dNx(2,:);
      B(3,1:2:16-1)=dNx(2,:);
      B(3,2:2:16)  =dNx(1,:);

      e=B*ed(ied,:)';
      if colD>3 ee([1 2 4])=e; else ee=e; end
      
      et=[et; ee'];
      es=[es; (D*ee)'];
    end
    
    ie=ie+incie;
  end
  
else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
