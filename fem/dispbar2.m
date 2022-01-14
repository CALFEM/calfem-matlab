function [sfac]=dispbar2(ex,ey,ed,plotpar,sfac)
%dispbar2(ex,ey,ed,plotpar,sfac)
%[sfac]=dispbar2(ex,ey,ed)
%[sfac]=dispbar2(ex,ey,ed,plotpar)
%
%------------------------------------------------------------------------
% PURPOSE: 
%  Draw the displacements for a two dimensional bar element.
%  	
% INPUT:
%	ex = [ x1 x2 ]
%	ey = [ y1 y2 ]	element node coordinates.
%
%   ed:     element displacement matrix
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

 c=size(ed);
 
 if (c(1)~=a(1)|c(2)~=4)
    disp('??? Check size of displacement input arguments!')
    return 
 end
 
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
   dx=ex(2)-ex(1);
   dy=ey(2)-ey(1);
   Length=sqrt(dx*dx+dy*dy);
   sfac=(0.1*Length)/max(max(abs(ed)));
 end
        
 x=(ex+sfac*ed(:,[1 3]))'; 
 y=(ey+sfac*ed(:,[2 4]))';
 xc=x;
 yc=y;
% ************* plot commands *******************
 hold on
 axis equal
 % Plots diagram.
 plot(xc,yc,s1)
 if plotpar(3)~=0
    plot(x,y,s2);
 end
 hold off; 
end

