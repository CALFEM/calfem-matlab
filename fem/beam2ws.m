 function es=beam2ws(ex,ey,ep,ed,eq)
% es=beam2ws(ex,ey,ep,ed,eq)
%-------------------------------------------------------------------------
%    PURPOSE
%      Compute section forces in a two dimensional beam element
%      on elastic foundation. 
% 
%    INPUT:  ex = [x1 x2]
%            ey = [y1 y2]            element node coordinates
%
%            ep = [E A I ka kt]      element properties,
%                                     E: Young's modulus
%                                     A: cross section area
%                                     I: moment of inertia
%                                     ka: axial foundation stiffness
%                                     kt: transversal found. stiffness
%
%            ed = [u1 ... u6]       element displacement vector
%
%            eq = [qx qy]           distributed loads, local directions
%          
%    OUTPUT: es = [N1 V1 M1 ;
%                  N2 V2 M2 ]       element forces, local directions 
% -------------------------------------------------------------------------
 
% LAST MODIFIED: K Persson  1996-04-26
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%--------------------------------------------------------------------------
  if length(ed(:,1)) > 1 
     disp('Only one row is allowed in the ed matrix !!!')
     return
  end
  
  b=[ ex(2)-ex(1); ey(2)-ey(1) ];
  L=sqrt(b'*b);  n=b/L;
  
  if nargin==4;  qx=0; qy=0; end
  if nargin==5;  qx=eq(1); qy=eq(2); end

  E=ep(1); A=ep(2); I=ep(3); ka=ep(4); kt=ep(5);

  K1=[ E*A/L   0            0      -E*A/L      0          0    ;
        0   12*E*I/L^3   6*E*I/L^2   0   -12*E*I/L^3  6*E*I/L^2;
        0   6*E*I/L^2    4*E*I/L     0   -6*E*I/L^2   2*E*I/L  ;
      -E*A/L   0            0       E*A/L      0          0    ;
        0   -12*E*I/L^3 -6*E*I/L^2   0   12*E*I/L^3  -6*E*I/L^2;
        0   6*E*I/L^2    2*E*I/L     0   -6*E*I/L^2   4*E*I/L ];

  K2=L/420*[140*ka   0       0      70*ka    0       0      ;
              0    156*kt   22*kt*L    0    54*kt  -13*kt*L ;
              0    22*kt*L  4*kt*L^2   0   13*kt*L -3*kt*L^2;   
             70*ka    0       0     140*ka    0       0     ;   
              0    54*kt    13*kt*L    0   156*kt  -22*kt*L ;                
              0   -13*kt*L -3*kt*L^2   0  -22*kt*L  4*kt*L^2];

  Kle=K1+K2;

  fle=L*[qx/2 qy/2 qy*L/12 qx/2 qy/2 -qy*L/12]';
  
  G=[n(1) n(2)  0    0    0   0;
    -n(2) n(1)  0    0    0   0;
      0    0    1    0    0   0;
      0    0    0   n(1) n(2) 0;
      0    0    0  -n(2) n(1) 0;
      0    0    0    0    0   1];
 
  P=(Kle*G*ed'-fle); 
  es=[-P(1) -P(2) -P(3) 
       P(4)  P(5)  P(6)];
%--------------------------------- end -------------------------------------
 
