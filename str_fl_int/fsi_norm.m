function  [XL,XR] = fsi_normc(XL,XR,laS,fiS,listS,laF,fiF,listF,H,fp)
%  [XL,XR] = fsi_normc(XL,XR,laS,fiS,listS,laF,fiF,listF,H,fp)
%-------------------------------------------------------------
% PURPOSE 
%   Mass normalise the left and right eigenvectors of 
%   the coupled problem (p-formulation). 
%
%   To be used between fsi_mod and fsi_egv!
%
%   
% INPUT:
%   XL :    left eigenvectors in reduced coord
%   XR :    right eigenvectors in reduced coord
%   laS :    column matrix, contains the eigenvalues, structure
%   fiS :    eigenvectors, structure
%   listS :  list of egenvectors, structure
%   laF :    column matrix, contains the eigenvalues, fluid
%   fiF :    eigenvectors, fluid
%   listF :  list of egenvectors, fluid  
%   H     :  coupling matrix
%   fp   :   fluid properties  [rho c]
%
% OUTPUT:
%   XL :    left eigenvectors in reduced coord mass normalized
%   XR :    right eigenvectors in reduced coord mass normalized
raa=fp(1);c=fp(2);
As=inv(diag(sqrt(raa*c^2*laS(listS))));
Af=inv(diag(sqrt(raa*c^2*laF(listF))));

Hstj=Af*fiF(:,listF)'*raa*c^2*H'*fiS(:,listS)*As;
Mstj=[As zeros(length(listS),length(listF));Hstj Af];
Mnorm=XL'*Mstj*XR;
for i=1:min(size(Mnorm))
	if real(Mnorm(i,i))<0
                XL(:,i)=-XL(:,i);
                XR(:,i)=XR(:,i);
	else
                XL(:,i)=XL(:,i);
                XR(:,i)=XR(:,i);
	end
end
Mnorm=XL'*Mstj*XR;
for i=1:min(size(Mnorm))
  d=inv((sqrt(Mnorm(i,i))));
  XL(:,i)=d'*XL(:,i);
  XR(:,i)=d*XR(:,i);
end
