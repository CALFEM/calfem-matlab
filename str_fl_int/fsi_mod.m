function  [Amodc] = fsipr(laS,fiS,listS,laF,fiF,listF,H,fp,string)
% [Amodc] = fsipr(laS,fiS,listS,laF,fiF,listF,H,fp,string)
%-------------------------------------------------------------
% PURPOSE 
%   Generates a coupled system of equations 
%   according to Eq. 32 or Eq. 32
%   
% INPUT:
%   laS :    column matrix, contains the eigenvalues, structure
%   fiS :    eigenvectors, structure
%   listS :  list of egenvectors, structure
%   laF :    column matrix, contains the eigenvalues, fluid
%   fiF :    eigenvectors, fluid
%   listF :  list of egenvectors, fluid  
%   H     :  coupling matrix
%   fp   :   fluid properties  [rho c]
%   string:  'right' or 'left'
%
% OUTPUT:
%   Amodc   coupled system matrix according to Eq. 31 or Eq. 32
% --------------------------------------------------
%
% Literature:
% Sandberg, G.: 
%	A new strategy for solving fluid-structure problems, 
%	International Journal for Numerical Methods in Engineering,
%	38, (1995), 357--370.
% --------------------------------------------------

% REFERENCES
%    G Sandberg  1996-08-07
% Copyright (c) 1991-96 by Division of Structural Mechanics and
%                          Department of Solid Mechanics.
%                          Lund Institute of Technology
%-------------------------------------------------------------
% ---- fluid properties ----------------------

rho=fp(1); c=fp(2); rhoc2=rho*c*c;

% ---- right system --------------------------

if string(1) == 'r'

% ---- structure -----------------------------

AmodS = diag(laS(listS));

% ---- coupling ------------------------------

Hmod=fiS(:,listS)'*H*fiF(:,listF);
Hmodc=-sqrt(rhoc2)*sqrt(AmodS)*Hmod;

% ---- fluid ---------------------------------

AmodF = diag(laF(listF)) + rhoc2*Hmod'*Hmod;

% ---- total system --------------------------

Amodc=[AmodS  Hmodc;
       Hmodc'  AmodF];
end

% ---- left system ---------------------------

if string(1) == 'l'

% ---- fluid ---------------------------------

AmodF = diag(laF(listF));

% ---- coupling ------------------------------

Hmod=fiS(:,listS)'*H*fiF(:,listF);
Hmodc=-sqrt(rhoc2)*Hmod*sqrt(AmodF);

% ---- structure -----------------------------

AmodS = diag(laS(listS)) + rhoc2*Hmod*Hmod';

% ---- total system --------------------------

Amodc=[AmodS  Hmodc;
       Hmodc'  AmodF];

end

%--------------------------end--------------------------------
