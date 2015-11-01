function [he]=cp4s4f(ex,ey,ez,ep)
% [he]=cp4s4f(ex,ey,ez,ep)
%-----------------------------------------------------------
% PURPOSE
%  Compute element coupling matrix between a 8 node 
%  linear acoustic element and a 8 node linear solid element. 
%
% INPUT:  ex = [x1 x2 x3 x4]   element coordinates
%         ey = [y1 y2 y3 y4]
%         ez = [z1 z2 z3 z4]
%                             
%         ep = [ir]   integration rule
%
% OUTPUT: he :  element coupling matrix (12 x 4)
%-------------------------------------------------------------

% LAST MODIFIED: P Davidsson    1998-10-20
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  ir=ep(1); ngp=ir*ir;

% this section for rotating the surface to a local x-y-plane 

  p2=[ ex(2)-ex(1); ey(2)-ey(1); ez(2)-ez(1) ];
  p3=[ ex(3)-ex(1); ey(3)-ey(1); ez(3)-ez(1) ];
  p4=[ ex(4)-ex(1); ey(4)-ey(1); ez(4)-ez(1) ];

  L=sqrt(p2'*p2);  n1=p2/L;
  L=sqrt(p4'*p4);  v2=p4/L;

  n3=[n1(2)*v2(3)-v2(2)*n1(3);
      n1(3)*v2(1)-v2(3)*n1(1);
      n1(1)*v2(2)-v2(1)*n1(2)];

  n3=n3/sqrt(n3'*n3);

  n2=[n3(2)*n1(3)-n1(2)*n3(3);
      n3(3)*n1(1)-n1(3)*n3(1);
      n3(1)*n1(2)-n1(1)*n3(2)];

  An= [n1,n2,n3];

  lp2=An'*p2;
  lp3=An'*p3;
  lp4=An'*p4;

  ec=[zeros(3,1),lp2,lp3,lp4];
  lex=ec(1,:);
  ley=ec(2,:);

  G=[n3', zeros(1,9); 
     zeros(1,3), n3', zeros(1,6);
     zeros(1,6), n3', zeros(1,3);
     zeros(1,9), n3']';

% this section for rotating ... ends here
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

  he1=zeros(4,4);  
  JT=dNr*[lex;ley]';

  for i=1:ngp
    indx=[ 2*i-1; 2*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminanten lika med noll!')
    end
    JTinv=inv(JT(indx,:));
    he1=he1+N(i,:)'*N(i,:)*detJ*wp(i);
  end

 he=G*he1;
%--------------------------end--------------------------------
