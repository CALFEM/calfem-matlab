  function es=beam2gs(ex,ey,ep,ed,N,eq)
% es=beam2gs(ex,ey,ep,ed,N,eq)
%-------------------------------------------------------------
%    PURPOSE
%      Calculate section forces in a two dimensional nonlinear
%      beam element.
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]          element node coordinates
%
%            ep = [E A I]          element properties,
%                                   E: Young's modulus
%                                   A: cross section area
%                                   I: moment of inertia
%
%            ed = [u1 ... u6]       element displacement vector
%
%            N			    axial force
% 
%            eq = [qy]              distributed transverse load
%          
%    OUTPUT: es = [N1 V1 M1 ;
%                  N2 V2 M2 ]       element forces, local directions 
%-------------------------------------------------------------

% LAST MODIFIED: K-G Olsson   1995-10-21
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
   b=[ ex(2)-ex(1); ey(2)-ey(1) ];
   L=sqrt(b'*b);  n=b/L;
% 
   if nargin==5; eq=0; end
%
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

   u=ed';
   P=(Kle*G*u-fle); 
   es=[-P(1) -P(2) -P(3) 
        P(4)  P(5)  P(6)];
%--------------------------end--------------------------------
