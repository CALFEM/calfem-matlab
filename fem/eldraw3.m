function eldraw3(ex,ey,ez,plotpar,elnum)
%eldraw3(ex,ey,ez,plotpar,elnum)
%eldraw3(ex,ey,ez,plotpar)
%eldraw3(ex,ey,ez)
%-------------------------------------------------------------
% PURPOSE 
%   Draw the undeformed 3D mesh for a number of elements of 
%   the same type. Supported elements are:
% 
%   1) -> bar element              2) -> beam el. 
%    
%   next:[ 3) -> tetraheder 3 node el.    4) -> brick 8 node el.]
%
% INPUT 
%    ex,ey,ez:.......... nen:   number of element nodes
%                        nel:   number of elements   
%    plotpar=[ linetype, linecolor, nodemark]
% 
%             linetype=1 -> solid    linecolor=1 -> black
%                      2 -> dashed             2 -> blue
%                      3 -> dotted             3 -> magenta
%                                              4 -> red
%
%             nodemark=1 -> circle       
%                      2 -> star               
%                      0 -> no mark 
%
%    elnum=edof(:,1) ;i.e. the first column in the topology matrix            
%               
%    Rem. Default is solid white lines with circles at nodes.
%         
%-------------------------------------------------------------

% LAST MODIFIED: A Olsson 2004-09-03
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==3)|(nargin==4)|(nargin==5))
    disp('??? Wrong number of input arguments!')
    return
 end 
 
 a=size(ex); b=size(ey); c=size(ez);
 
 if ((a-b)==[0 0])&((b-c)==[0 0])
    nel=a(1);nen=a(2);
 else
    disp('??? Check size of coordinate input arguments!')
    return
 end
 if nargin==3; 
      plotpar=[1 1 1]; 
 end
 
 [s1,s2]=pltstyle(plotpar);
 
% ************************************************** 
% ************* plot coordinates *******************
% ************************************************** 
 
 x0=sum(ex')/nen; y0=sum(ey')/nen; z0=sum(ez')/nen;
  
% ********** Bar or Beam elements *************
 if nen==2
    x=ex'; y=ey'; z=ez';  
    xc=x; yc=y; zc=z;
%**********************************************************          
 else
    disp('!!! Sorry, this element is currently not supported!')
    return      
 end
  
%*************************************************
% ************** plot commands *******************
%*************************************************
 axis('equal')
 hold on  
 view(3)
 plot3(xc,yc,zc,s1)   
 plot3(x,y,z,s2)
 if nargin==5
    for i=1:nel
       text(x0(i),y0(i),z0(i),int2str(elnum(i)))
    end 
 end
 xlabel('x'); ylabel('y'); zlabel('z');
 hold off 
%--------------------------end-------------------------------- 
