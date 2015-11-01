function [Ke,fe]=platre(ex,ey,ep,D,eq)
% Ke=platre(ex,ey,ep,D)
% [Ke,fe]=platre(ex,ey,ep,D,eq)
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a rectangular plate element.
%  NOTE! Element sides must be parallel to the coordinate axis.
%
% INPUT:  ex = [x1 x2 x3 x4]      element coordinates
%         ey = [y1 y2 y3 y4]
%
%         ep=[t]                  thicknes
%
%         D                       constitutive matrix for
%                                 plane stress         
%      
%         eq=[qz]                 load/unit area
%
%  OUTPUT: Ke :  element stiffness matrix (12 x 12)
%          fe : equivalent nodal forces (12 x 1)
%-------------------------------------------------------------

% LAST MODIFIED: K Persson   1995-08-23
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
Lx=ex(3)-ex(1); Ly=ey(3)-ey(1); t=ep(1,1);
%
D=t^3/12*D;
%
A1=Ly/(Lx^3); A2=Lx/(Ly^3); A3=1/Lx/Ly;
A4=Ly/(Lx^2); A5=Lx/(Ly^2); A6=1/Lx;
A7=1/Ly;      A8=Ly/Lx;     A9=Lx/Ly;
%
C1= 4*A1*D(1,1)+4*A2*D(2,2)+2*A3*D(1,2)+5.6*A3*D(3,3);
C2=-4*A1*D(1,1)+2*A2*D(2,2)-2*A3*D(1,2)-5.6*A3*D(3,3);
C3= 2*A1*D(1,1)-4*A2*D(2,2)-2*A3*D(1,2)-5.6*A3*D(3,3);
C4=-2*A1*D(1,1)-2*A2*D(2,2)+2*A3*D(1,2)+5.6*A3*D(3,3);
C5=2*A5*D(2,2)+A6*D(1,2)+0.4*A6*D(3,3);
C6=2*A4*D(1,1)+A7*D(1,2)+0.4*A7*D(3,3);

C7=2*A5*D(2,2)+0.4*A6*D(3,3);
C8=2*A4*D(1,1)+0.4*A7*D(3,3);
C9= A5*D(2,2)-A6*D(1,2)-0.4*A6*D(3,3);
C10=A4*D(1,1)-A7*D(1,2)-0.4*A7*D(3,3);
C11=A5*D(2,2)-0.4*A6*D(3,3);
C12=A4*D(1,1)-0.4*A7*D(3,3);

C13=4/3*A9*D(2,2)+8/15*A8*D(3,3);
C14=4/3*A8*D(1,1)+8/15*A9*D(3,3);
C15=2/3*A9*D(2,2)-8/15*A8*D(3,3);
C16=2/3*A8*D(1,1)-8/15*A9*D(3,3);
C17=2/3*A9*D(2,2)-2/15*A8*D(3,3);
C18=2/3*A8*D(1,1)-2/15*A9*D(3,3);
C19=1/3*A9*D(2,2)+2/15*A8*D(3,3);
C20=1/3*A8*D(1,1)+2/15*A9*D(3,3);
C21=D(1,2);
%
Keq=zeros(12,12);
Keq(1,1:12)=[C1 C5 -C6 C2 C9 -C8 C4 C11 -C12 C3 C7 -C10];
Keq(2,2:12)=[C13 -C21 C9 C15 0 -C11 C19 0 -C7 C17 0];
Keq(3,3:12)=[C14 C8 0 C18 C12 0 C20 -C10 0 C16];
Keq(4,4:12)=[C1 C5 C6 C3 C7 C10 C4 C11 C12];
Keq(5,5:12)=[C13 C21 -C7 C17 0 -C11 C19 0];
Keq(6,6:12)=[C14 C10 0 C16 -C12 0 C20];
Keq(7,7:12)=[C1 -C5 C6 C2 -C9 C8];
Keq(8,8:12)=[C13 -C21 -C9 C15 0];
Keq(9,9:12)=[C14 -C8 0 C18];
Keq(10,10:12)=[C1 -C5 -C6];
Keq(11,11:12)=[C13 C21];
Keq(12,12)=[C14];
Keq=Keq'+Keq-diag(diag(Keq));
%
if nargin==5
R1=eq*Lx*Ly/4;
R2=eq*Lx*Ly^2/24;
R3=eq*Ly*Lx^2/24;
%
feq(:,1)=[R1 R2 -R3 R1 R2 R3 R1 -R2 R3 R1 -R2 -R3]';
fe=feq;
end
Ke=Keq;  
%--------------------------end--------------------------------
