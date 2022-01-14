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

% LAST MODIFIED: Kent Persson   2009-11-30
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%-------------------------------------------------------------
%
Lx=abs(ex(3)-ex(1)); Ly=abs(ey(3)-ey(1));t=ep(1,1);
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
ed1=ed(:,1); ed2=ed(:,2); ed3=ed(:,3);
ed4=ed(:,4); ed5=ed(:,5); ed6=ed(:,6);
ed7=ed(:,7); ed8=ed(:,8); ed9=ed(:,9);
ed10=ed(:,10); ed11=ed(:,11); ed12=ed(:,12);
%
if ((ex(3)-ex(1))<0) & (ey(3)-ey(1)>0)
ed1=ed(:,10); ed2=ed(:,11); ed3=ed(:,12);
ed4=ed(:,1); ed5=ed(:,2); ed6=ed(:,3);
ed7=ed(:,4); ed8=ed(:,5); ed9=ed(:,6);
ed10=ed(:,7); ed11=ed(:,8); ed12=ed(:,9);
end    
%
if ((ex(3)-ex(1))<0) & (ey(3)-ey(1)<0)
ed1=ed(:,7); ed2=ed(:,8); ed3=ed(:,9);
ed4=ed(:,10); ed5=ed(:,11); ed6=ed(:,12);
ed7=ed(:,1); ed8=ed(:,2); ed9=ed(:,3);
ed10=ed(:,4); ed11=ed(:,5); ed12=ed(:,6);
end    
%
if ((ex(3)-ex(1))>0) & (ey(3)-ey(1)<0)
ed1=ed(:,4); ed2=ed(:,5); ed3=ed(:,6);
ed4=ed(:,7); ed5=ed(:,8); ed6=ed(:,9);
ed7=ed(:,10); ed8=ed(:,11); ed9=ed(:,12);
ed10=ed(:,1); ed11=ed(:,2); ed12=ed(:,3);
end    
%
mx=A3*(-ed2-ed5+ed8+ed11)+A2*(ed3-ed6-ed9+ed12);
my=A1*(-ed2-ed5+ed8+ed11)+A4*(ed3-ed6-ed9+ed12);
mxy=A6*(ed2-ed5-ed8+ed11)+A5*(-ed3-ed6+ed9+ed12)+A7*(ed1-ed4+ed7-ed10);

vx=B5*(-ed2+ed5-ed8+ed11)+B4*(ed3+ed6+ed9+ed12)+B2*(-ed1+ed4+ed7-ed10);
vy=B3*(ed2+ed5+ed8+ed11)+B5*(ed3-ed6+ed9-ed12)+B1*(-ed1-ed4+ed7+ed10);
     
m1=0.5.*(mx+my)+sqrt(0.25.*(mx-my).^2+mxy.^2);
m2=0.5.*(mx+my)-sqrt(0.25.*(mx-my).^2+mxy.^2);
alfa=0.5*180/pi*atan2(mxy,(mx-my)/2);
         
es=[mx my mxy vx vy];

et=-inv(D)*[mx;my;mxy];
et=et';
%--------------------------end--------------------------------