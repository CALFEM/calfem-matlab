 function es=beam2ds(ex,ey,ep,ed,ev,ea)
% es=beam2ds(ex,ey,ep,ed,ev,ea)
%-------------------------------------------------------------
%    PURPOSE
%     Calculate the element forces for a number of identical 
%     (nie) 2D Bernoulli beam elements in dynamic analysis. 
%
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]         element node coordinates
%
%            ep = [E A I m (a b)]
%                                 E: Young's modulus
%                                 A: cross section area
%                                 I: moment of inertia
%                                 m: mass per unit length
%                               a,b: damping coefficients
%                                    Ce=a*Me+b*Ke
%
%             ed :  element displacement matrix 
%             ev :  element velocity matrix 
%             ea :  element acceleration matrix 
%          
%    OUTPUT: es : element forces in local directions,
%               = [-N1 -V1 -M1 N2 V2 M2;
%                  .......      ......] ; dim(es)= nie x 6
%-------------------------------------------------------------
 
% LAST MODIFIED: K Persson    1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;
% 
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
   [nie,ned]=size(ed);
   for i=1:nie
        d=ed(i,:)';
        v=ev(i,:)';
        a=ea(i,:)';
        es(i,:)=(Kle*G*d+Cle*G*v+Mle*G*a)';
   end
%--------------------------end--------------------------------
