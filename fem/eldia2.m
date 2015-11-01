function [sfac]=eldia2(ex,ey,es,plotpar,sfac,eci)
%eldia2(ex,ey,es,plotpar,sfac)
%eldia2(ex,ey,es,plotpar,sfac,eci)
%[sfac]=eldia2(ex,ey,es)
%[sfac]=eldia2(ex,ey,es,plotpar)
%
%-------------------------------------------------------------
% PURPOSE: 
%  Draw the section force diagrams of a two dimensional beam element.
%  	
% INPUT:
%	ex = [ x1 x2 ]
%	ey = [ y1 y2 ]	element node coordinates.
%
%	es = [  X1;
%		X2;
%		... ] 	vector containing the section force
%			in Nbr evaluation points along the beam.
%	
%    plotpar=[ linecolor, elementcolor] 
%
%             linecolor=1 -> black    elementcolor=1 -> black
%                       2 -> blue                    2 -> blue
%                       3 -> magenta                 3 -> magenta
%                       4 -> red                     4 -> red
%
%	sfac = [scalar]	scale factor.
%
%	eci = [	x1;
%		x2;
%		... ]	local x-coordinates of the evaluation points (Nbr).
%       If not given, the evaluation points are assumed to be uniformly
%       distributed.
%	
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2004-10-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------

 if ~((nargin==3)|(nargin==4)|(nargin==5)|(nargin==6))
    disp('??? Wrong number of input arguments!')
    return
 end
 a=size(ex); b=size(ey);
 
 if ~((a-b)==[0 0])
    disp('??? Check size of coordinate input arguments!')
    return
 end

 c=size(es);
 Nbr=c(1);
 b=[ex(2)-ex(1);ey(2)-ey(1)];
 Length=sqrt(b'*b);
 n=b/Length;
 if nargin==3
    plotpar=[2 1]
 end
 if((nargin==3)|(nargin==4))
    sfac=(0.2*Length)/max(abs(es));
 end
 if((nargin==3)|(nargin==4)|(nargin==5))
    eci=[0:Length/(Nbr-1):Length]'; 
 end   
 if (size(plotpar)==[1 2])
    if plotpar(1)==1 ; s1='-k';
    elseif plotpar(1)==2 ; s1='-b';
    elseif plotpar(1)==3 ; s1='-m';
    elseif plotpar(1)==4 ; s1='-r';
    else disp('??? Error in variable plotpar(1)!');
       return;
    end
    if plotpar(2)==1 ; s2='-k';
    elseif plotpar(2)==2 ; s2='-b';
    elseif plotpar(2)==3 ; s2='-m';
    elseif plotpar(2)==4 ; s2='-r';
    else disp('??? Error in variable plotpar(2)!');
       return;
    end
 else
    disp('??? Check size of "plotpar" input argument!')
    return
 end
 a=size(eci); 
 if ~(c(1)==a(1)|(c(2)==1))
    disp('??? Check size of "es" or "eci" input arguments!')
    return 
 end 
       
 es=es*sfac;
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
	A(i,1)=A(i,1)+es(i)*n(2);
	A(i,2)=A(i,2)-es(i)*n(1);
 end

% ************* plot commands *******************
 hold on
 axis equal
 % Plots diagram.
 plot(A(:,1),A(:,2),s1);
 plot([ex(1) A(1,1)],[ey(1) A(1,2)],s1);
 plot([ex(2) A(Nbr,1)],[ey(2) A(Nbr,2)],s1);
% Plots stripes in diagram.
 for i=1:Nbr
	plot([B(i,1) A(i,1)],[B(i,2) A(i,2)],s1);
 end
% Plots element.	
 a=plot(ex',ey',s2);
 set(a,'LineWidth',[2]);
 hold off; 
end

