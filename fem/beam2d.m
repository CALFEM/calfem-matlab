function [Ke,Me,Ce]=beam2d(ex,ey,ep)
% [Ke,Me]=beam2d(ex,ey,ep)
% [Ke,Me,Ce]=beam2d(ex,ey,ep)
%-------------------------------------------------------------
%   PURPOSE
%    Calculate the stiffness matrix Ke, the mass matrix Me 
%    and the damping matrix Ce for a 2D elastic Bernoulli
%    beam element. 
%
%   INPUT:   ex = [x1 x2]
%            ey = [y1 y2]         element node coordinates 
%
%            ep = [E A I m (a b)]  
%                                 E: Young's modulus
%                                 A: cross section area
%                                 I: moment of inertia
%                                 m: mass per unit length
%                               a,b: damping coefficients,
%                                    Ce=aMe+bKe 
% 
%   OUTPUT:  Ke : element stiffness matrix (6 x 6)
%            Me : element mass matrix 
%            Ce : element damping matrix, optional
%-------------------------------------------------------------
 
% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;
   
   E=ep(1);    A=ep(2);    I=ep(3);     m=ep(4);
   a=0 ; b=0 ;
   if length(ep)==6 ; a=ep(5) ; b=ep(6) ; end
%
   Kle=[E*A/L   0            0      -E*A/L      0          0 ;
          0   12*E*I/L^3   6*E*I/L^2  0   -12*E*I/L^3  6*E*I/L^2;
          0   6*E*I/L^2    4*E*I/L    0   -6*E*I/L^2   2*E*I/L;
        -E*A/L  0            0       E*A/L      0          0 ;
          0   -12*E*I/L^3 -6*E*I/L^2  0   12*E*I/L^3  -6*E*I/L^2;
          0   6*E*I/L^2    2*E*I/L    0   -6*E*I/L^2   4*E*I/L];
%
   Mle=m*L/420*[140   0     0    70   0    0    ;
                 0   156   22*L   0   54  -13*L ;
                 0   22*L  4*L^2  0  13*L -3*L^2 ;
                70    0     0   140   0     0   ;
                 0   54    13*L   0  156  -22*L ;
                 0  -13*L -3*L^2  0 -22*L  4*L^2];
%
   Cle=a*Mle+b*Kle;
%
   G=[n(1) n(2)  0    0    0   0;
     -n(2) n(1)  0    0    0   0;
       0    0    1    0    0   0;
       0    0    0   n(1) n(2) 0;
       0    0    0  -n(2) n(1) 0;
       0    0    0    0    0   1];
%
   Ke=G'*Kle*G;    Me=G'*Mle*G;    Ce=G'*Cle*G;
%--------------------------end--------------------------------
 
