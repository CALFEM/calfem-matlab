function pltscalb2(sfac,magnitude,plotpar)
%pltscalb2(sfac,magnitude)
%pltscalb2(sfac,magnitude,plotpar)
%
%-------------------------------------------------------------
% PURPOSE: 
%  Draw a scale bar
%  	
% INPUT:
%
%	sfac = [scalar]	scale factor.
%	
%	magnitude = [Ref x y]	The scale bar has a length equivalent   
%				to Ref and starts at coordinates (x,y).
%				If no coordinates are given the starting
%				point will be (0,-0.5).
%
%    plotpar=[linecolor]
% 
%             linecolor=1 -> black
%             2 -> blue
%             3 -> magenta
%             4 -> red
%
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2004-09-13
% Copyright (c)  Division of Structural Mechanics
%                Lund Institute of Technology
%-------------------------------------------------------------

if ~((nargin==2)|(nargin==3))
    disp('??? Wrong number of input arguments!')
    return
end
fac=size(sfac);
mag=size(magnitude);
if ~((fac==[1 1]))
   disp('??? Check size of coordinate input arguments!')
   return
end
if ~((mag==[1 1])|(mag==[1 3]))
   disp('??? Check size of coordinate input arguments!')
   return
end

hold on
% Creating a scale bar.
N=magnitude(1);
L=N*sfac;
mag=size(magnitude);
if mag(2)==1
	x=0;
	y=-0.5;
else
	x=magnitude(2);
	y=magnitude(3);
end;

if ((nargin==2))
    plotpar=2;
end
if plotpar==1 ; s1=['k'];
 elseif plotpar==2 ; s1=['b'];
 elseif plotpar==3 ; s1=['m'];
 elseif plotpar==4 ; s1=['r'];
 else disp('??? Error in variable plotpar!');
      return;
 end

plot([x (x+L)],[y y],s1);
plot([x x],[(y-L/20) (y+L/20)],s1);
%line([(x+L/2) (x+L/2)],[(y-L/20) (y+L/20)]);
plot([(x+L) (x+L)],[(y-L/20) (y+L/20)],s1);
text((N*1.1)*sfac+x,y,sprintf('%0.5g ', N));

hold off;

