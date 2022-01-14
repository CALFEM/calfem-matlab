function [sfac]=dispbeam2(ex,ey,edi,plotpar,sfac)
%dispbeam2(ex,ey,edi,plotpar,sfac)
%[sfac]=dispbeam2(ex,ey,edi)
%[sfac]=dispbeam2(ex,ey,edi,plotpar)
%
%------------------------------------------------------------------------
% PURPOSE: 
%  Draw the displacement diagram for a two dimensional beam element.
%  	
% INPUT:
%	ex = [ x1 x2 ]
%	ey = [ y1 y2 ]	element node coordinates.
%
%	edi = [ u1 v1;
%		u2 v2;
%		... ] 	matrix containing the displacements
%			in Nbr evaluation points along the beam.
%	
%    plotpar=[linetype, linecolour, nodemark] 
%
%             linetype=1 -> solid   linecolour=1 -> black
%                      2 -> dashed             2 -> blue
%                      3 -> dotted             3 -> magenta
%                                              4 -> red
%             nodemark=0 -> no mark 
%                      1 -> circle       
%                      2 -> star              
%                      3 -> point 
%
%    sfac = [scalar] scale factor for displacements. 
%            
%    Rem. Default if sfac and plotpar is left out is auto magnification 
%         and dashed black lines with circles at nodes -> plotpar=[1 1 1] 
%------------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2015-11-18
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%------------------------------------------------------------------------

 if ~((nargin==3)|(nargin==4)|(nargin==5))
    disp('??? Wrong number of input arguments!')
    return
 end
 a=size(ex); b=size(ey);
 
 if ~((a-b)==[0 0])
    disp('??? Check size of coordinate input arguments!')
    return
 end

 c=size(edi);
 if (c(2)~=2)
    disp('??? Check size of displacement input argument!')
    return
 end
 Nbr=c(1);
 b=[ex(2)-ex(1);ey(2)-ey(1)];
 Length=sqrt(b'*b);
 n=b/Length;
 
 if nargin==3
    plotpar=[1 1 1]
 end
 if (size(plotpar)==[1 3])
     if plotpar(1)==1 ; lt='-';
     elseif plotpar(1)==2 ; lt='--';
     elseif plotpar(1)==3 ; lt=':';
     else disp('??? Error in variable plotpar(1)!');
          return;
     end

     if plotpar(2)==1 ; lc='k';
     elseif plotpar(2)==2 ; lc='b';
     elseif plotpar(2)==3 ; lc='m';
     elseif plotpar(2)==4 ; lc='r';
     else disp('??? Error in variable plotpar(2)!');
          return;
     end

     if plotpar(3)==0 ; nm='';
     elseif plotpar(3)==1 ; nm='o';
     elseif plotpar(3)==2 ; nm='*';
     elseif plotpar(3)==3 ; nm='.';    
     else disp('??? Error in variable plotpar(3)!');
          return;
     end

     s1=[lt lc];
     s2=[lc nm];

 else
   disp('??? Check size of "plotpar" input argument!')
   return
 end
 
 if((nargin==3)|(nargin==4))
    sfac=(0.1*Length)/max(max(abs(edi)));
 end
 eci=[0:Length/(Nbr-1):Length]'; 
        
 edi=edi*sfac;
% From local x-coordinates to global coordinates of the beam element.
 A=zeros(Nbr,2);
 A(1,1)=ex(1);
 A(1,2)=ey(1);
 for i=2:Nbr
	A(i,1)=A(1,1)+eci(i)*n(1);
	A(i,2)=A(1,2)+eci(i)*n(2);
 end

 B=A;
 for i=1:Nbr
	A(i,1)=A(i,1)+edi(i,1)*n(1)-edi(i,2)*n(2);
	A(i,2)=A(i,2)+edi(i,1)*n(2)+edi(i,2)*n(1);
 end
 A1=[A(1,1) A(Nbr,1)];
 A2=[A(1,2) A(Nbr,2)];
% ************* plot commands *******************
 hold on
 axis equal
 % Plots diagram.
 a=plot(A(:,1),A(:,2),s1);
 if plotpar(3)~=0
 a=plot(A1,A2,s2);
 end
% plot([ex(1) A(1,1)],[ey(1) A(1,2)],s1);
% plot([ex(2) A(Nbr,1)],[ey(2) A(Nbr,2)],s1);
% Plots stripes in diagram.
% for i=1:Nbr
%	plot([B(i,1) A(i,1)],[B(i,2) A(i,2)],s1);
% end
% Plots element.	
% plot(ex',ey',s2);
 
 hold off; 
end

