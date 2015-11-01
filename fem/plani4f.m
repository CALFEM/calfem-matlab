 function [ef]=plani4f(ex,ey,ep,es)
% ef=plani4f(ex,ey,ep,es)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the element force vector corresponding to the 
%  stresses in a 4 node isoparametric element.
%
% INPUT:  ex = [x1 x2 x3 x4]       element coordinates
%         ey = [y1 y2 y3 y4]
%             
%         ep = [ptype t ir]        ptype: analysis type
%                                  ir: integration rule
%                                  t : thickness
%
%         es = [ sigx sigy [sigz] tauxy  element stress matrix
%                  ......              ]  one row for each integration point
%
% OUTPUT: fe = [f1 f2 ...f8]';     internal force vector
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa  1995-10-25
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  ptype=ep(1); t=ep(2); ir=ep(3); ngp=ir*ir;
   
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

%--------- plane stress --------------------------------------
if ptype==1
   
  [rowes,colD]=size(es);
  rowex=size(ex,1);

  if rowex==1 incie=0; else incie=1; end
  
  ef=[]; ir=0; ie=1;
  for ied=1:rowes/ngp
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    fint=zeros(8,1);
    for i=1:ngp
      ir=ir+1;
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

      if colD>3 stress=es(ir,[1 2 4]); else stress=es(ir,:); end
      fint=fint+B'*stress'*wp(i)*detJ*t;
    end
    
   ef=[ef; fint'];
   ie=ie+incie;
  end
  
%--------- plane strain --------------------------------------
elseif ptype==2
   
  [rowes,colD]=size(es);
  rowex=size(ex,1);

  if rowex==1 incie=0; else incie=1; end
  
  ef=[]; ir=0; ie=1;
  res=size(es,1);
  for ied=1:rowes/ngp
    JT=dNr*[ex(ie,:);ey(ie,:)]';

    fint=zeros(8,1);
    for i=1:ngp
      ir=ir+1;
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

      if colD>3 stress=es(ir,[1 2 4]); else stress=es(ir,:); end
      fint=fint+B'*stress'*wp(i)*detJ*t;
    end
    
    ef=[ef; fint'];
    ie=ie+incie;
  end

else
   error('Error ! Check first argument, ptype=1 or 2 allowed')
   return
end
%--------------------------end--------------------------------
