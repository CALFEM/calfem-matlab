function [es,et]=platrs(ex,ey,ep,D,ed)
% [es,et]=platrs(ex,ey,ep,D,ed)
%-------------------------------------------------------------
% PURPOSE
%  Calculate element normal and shear force for a
%  rectangular plate element.
%
% INPUT: ex = [x1 x2 x3 x4]     element coordinates
%        ey = [y1 y2 y3 y4]
%
%        ep = [t]               thickness
%
%        D                      constitutive matrix for
%                               plane stress
%
%        ed = [u1 u2......u12;  element displacement vector
%             .............   ] one row for each element
%
% OUTPUT: es = [ Mxx Myy Mxy Vxz Vyz;   element force matrix
%                    ......          ]  one row for each element
%         et = [kxx,kyy,kxy]       curvature in global coordinates
%-------------------------------------------------------------

% LAST MODIFIED: Anders Olsson   1999-03-01
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
Lx=ex(3)-ex(1); Ly=ey(3)-ey(1);t=ep(1,1);
%
D=t^3/12*D;
%
A1=D(2,2)/2/Ly;
A2=D(1,1)/2/Lx;
A3=D(1,2)/2/Ly;
A4=D(1,2)/2/Lx;
A5=D(3,3)/2/Ly;
A6=D(3,3)/2/Lx;
A7=4*D(3,3)/Lx/Ly;

B1=6*D(2,2)/Ly/Ly/Ly;
B2=6*D(1,1)/Lx/Lx/Lx;
B3=-3*D(2,2)/Ly/Ly;
B4=3*D(1,1)/Lx/Lx;
B5=D(1,2)/Lx/Ly;
%
mx=A3*(-ed(:,2)-ed(:,5)+ed(:,8)+ed(:,11))+A2*(ed(:,3)-ed(:,6)-ed(:,9)+ed(:,12));
my=A1*(-ed(:,2)-ed(:,5)+ed(:,8)+ed(:,11))+A4*(ed(:,3)-ed(:,6)-ed(:,9)+ed(:,12));
mxy=A6*(ed(:,2)-ed(:,5)-ed(:,8)+ed(:,11))+A5*(-ed(:,3)-ed(:,6)+ed(:,9)+ed(:,12))...
     +A7*(ed(:,1)-ed(:,4)+ed(:,7)-ed(:,10));
     
m1=0.5.*(mx+my)+sqrt(0.25.*(mx-my).^2+mxy.^2);
m2=0.5.*(mx+my)-sqrt(0.25.*(mx-my).^2+mxy.^2);
alfa=0.5*180/pi*atan2(mxy,(mx-my)/2);

vx=B5*(-ed(:,2)+ed(:,5)-ed(:,8)+ed(:,11))+B4*(ed(:,3)+ed(:,6)+ed(:,9)+ed(:,12))...
     +B2*(-ed(:,1)+ed(:,4)+ed(:,7)-ed(:,10));
vy=B3*(ed(:,2)+ed(:,5)+ed(:,8)+ed(:,11))+B5*(ed(:,3)-ed(:,6)+ed(:,9)-ed(:,12))...
     +B1*(-ed(:,1)-ed(:,4)+ed(:,7)+ed(:,10));
         
         
es=[mx my mxy vx vy];

et=-inv(D)*[mx;my;mxy];
et=et';
%--------------------------end--------------------------------
