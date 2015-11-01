 function [Ke,fe]=beam2g(ex,ey,ep,N,eq)
% Ke=beam2g(ex,ey,ep,N)
% [Ke,fe]=beam2g(ex,ey,ep,N,eq)
%-------------------------------------------------------------
%    PURPOSE
%       Compute the element stiffness matrix for a two dimensional
%       beam element with respect to geometric nonlinearity.
%
%    INPUT: ex = [x1 x2]   element node coordinates
%           ey = [y1 y2]
%  
%           ep = [E A I]   element properties;
%                            E: Young's modulus
%                            A: Cross section area
%                            I: Moment of inertia
%
%           N:  axial force in the beam.
% 
%           eq: distributed transverse load
%              
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%
%            fe : element load vector (6 x 1)  
%-------------------------------------------------------------

% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

  if nargin==4; eq=0; end 
  
  if length(eq)>1
    disp('eq should be a scalar!!!')
    return 
  end

  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;

  E=ep(1); A=ep(2); I=ep(3); rho=-N*L^2/(pi^2*E*I);

  kL=pi*sqrt(abs(rho))+eps;

      if rho>0
         f1=(kL/2)/tan(kL/2);
         f2=(1/12)*kL^2/(1-f1);
         f3=f1/4+3*f2/4;
         f4=-f1/2+3*f2/2;
         f5=f1*f2;

         h=6*(2/kL^2-(1+cos(kL))/(kL*sin(kL)));

      elseif rho<0
         f1=(kL/2)/tanh(kL/2);
         f2=-(1/12)*kL^2/(1-f1);
         f3=f1/4+3*f2/4;
         f4=-f1/2+3*f2/2;
         f5=f1*f2;

         h=-6*(2/kL^2-(1+cosh(kL))/(kL*sinh(kL)));
      else
         f1=1;f2=1;f3=1;f4=1;f5=1;h=1;
      end

   Kle=[E*A/L  0            0        -E*A/L      0          0 ;
       0  12*E*I*f5/L^3   6*E*I*f2/L^2  0 -12*E*I*f5/L^3 6*E*I*f2/L^2;
       0  6*E*I*f2/L^2    4*E*I*f3/L    0  -6*E*I*f2/L^2   2*E*I*f4/L;
       -E*A/L  0            0         E*A/L      0          0 ;
       0  -12*E*I*f5/L^3 -6*E*I*f2/L^2  0  12*E*I*f5/L^3 -6*E*I*f2/L^2;
       0  6*E*I*f2/L^2    2*E*I*f4/L    0  -6*E*I*f2/L^2   4*E*I*f3/L];

   fle=eq*L*[0 1/2 L*h/12 0 1/2 -L*h/12]';


  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];

  Ke=G'*Kle*G;  fe=G'*fle;
%--------------------------end--------------------------------

 
