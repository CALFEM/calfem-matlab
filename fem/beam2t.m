function [Ke,fe]=beam2t(ex,ey,ep,eq)
%       function [Ke,fe]=beam2t(ex,ey,ep,eq)
%-------------------------------------------------------------
%    PURPOSE
%     Compute the stiffness matrix for a two dimensional elastic 
%     Timoshenko beam element. 
% 
%    INPUT:  ex = [x1 x2]      
%            ey = [y1 y2]      element node coordinates
%
%            ep = [E G A I ks ]     element properties,
%                               E: Young's modulus
%                               G: shear modulus
%                               A: Cross section area
%                               I: Moment of inertia
%                              ks: shear correction factor
%
%            eq = [qx qy]      distributed loads, local directions
% 
%    OUTPUT: Ke : element stiffness matrix (6 x 6)
%
%            fe : element load vector (6 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: E Serrano   1995-08-23 
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  if nargin==3;   eq=[0 0];    end
%
  b=[ex(2)-ex(1);ey(2)-ey(1)];
  L=sqrt(b'*b);  n=b/L;
  E=ep(1);
  Gm=ep(2);
  A=ep(3);
  I=ep(4);
  ks=ep(5);
%
  m=(12/L^2)*(E*I/(Gm*A*ks));
%
  Kle=E/(1+m)*[A*(1+m)/L    0         0       -A*(1+m)/L   0            0;
                  0      12*I/L^3  6*I/L^2        0      -12*I/L^3   6*I/L^2;
                  0      6*I/L^2  4*I*(1+m/4)/L   0      -6*I/L^2  2*I*(1-m/2)/L;
              -A*(1+m)/L    0         0        A*(1+m)/L    0           0;
                  0     -12*I/L^3  -6*I/L^2       0      12*I/L^3   -6*I/L^2;
                  0      6*I/L^2  2*I*(1-m/2)/L   0      -6*I/L^2   4*I*(1+m/4)/L];
%
  fle=L*[eq(1)/2  eq(2)/2  eq(2)*L/12  eq(1)/2  eq(2)/2  -eq(2)*L/12]';
%
  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
%
  Ke=G'*Kle*G;  fe=G'*fle;
%--------------------------- end -----------------------------
