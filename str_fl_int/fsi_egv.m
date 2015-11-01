function  [V] = fsi_egv(v,laS,fiS,listS,laF,fiF,listF,fp,string1,string2)
%  [V] = fsi_egv(v,laS,fiS,listS,laF,fiF,listF,fp,string1,string2)
%-------------------------------------------------------------
% PURPOSE 
%   calculates the right eigenvectors in the original FE
%   coordinate system, cf. Eqs. 33, 34, 35, 36 
%   
% INPUT:
%   v   :    eigenvectors in modal coordinates
%   laS :    column matrix, contains the eigenvalues, structure
%   fiS :    eigenvectors, structure
%   listS :  list of egenvectors, structure
%   laF :    column matrix, contains the eigenvalues, fluid
%   fiF :    eigenvectors, fluid
%   listF :  list of egenvectors, fluid  
%   fp   :   fluid properties  [rho c]
%
% OUTPUT:
%   V       eigenvectors in FE coord. syst.
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

rho=fp(1); c=fp(2); rhoc2=rho*c*c; c2rho=c*c/rho;

nolla=zeros(length(fiS),length(listF));
nollb=zeros(length(fiF),length(listS));

%-------------------------------------------------------------
if string2(1) == 'p'

  % ---- cf. Eq. 33 ---------------------------
  if string1(1) == 'r'

    af=fiS(:,listS)*diag((sqrt(rhoc2*laS(listS))).^(-1));
    V = [ af      nolla;
          nollb   fiF(:,listF) ]*v;

  end
  % -------------------------------------------

  % ---- cf. Eq. 34 ---------------------------
    if string1(1) == 'l'

    af=fiF(:,listF)*diag((sqrt(rhoc2*laF(listF))).^(-1));

    V = [ fiS(:,listS)   nolla;
          nollb          af   ]*v;

  end
  % -------------------------------------------
end
%-------------------------------------------------------------

%-------------------------------------------------------------
if string2(1) == 'd'

  % ---- cf. Eq. 35 ---------------------------
  if string1(1) == 'r'

    af=fiF(:,listF)*diag(sqrt(c2rho*(laF(listF)).^(-1)));

    V = [ fiS(:,listS)   nolla;
          nollb          af   ]*v;
    
  end
  % -------------------------------------------

  % ---- cf. Eq. 36 ---------------------------
  if string1(1) == 'l'

    af=fiS(:,listS)*diag(sqrt(c2rho*(laS(listS)).^(-1)));

    V = [ af      nolla;
          nollb   fiF(:,listF) ]*v;

  end
  % -------------------------------------------
end

%--------------------------end--------------------------------
