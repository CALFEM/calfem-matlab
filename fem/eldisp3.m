function [magnfac]=eldisp3(ex,ey,ez,ed,plotpar,magnfac)
%eldisp3(ex,ey,ez,ed,plotpar,magnfac)
%[magnfac]=eldisp3(ex,ey,ez,ed)
%-------------------------------------------------------------
% PURPOSE 
%   Draw the deformed 3D mesh for a number of elements of 
%   the same type. Supported elements are:
% 
%     1) -> bar element              2) -> beam el.
%  
%  to be supported:[ 3) -> tetraheder 4 node el.  4) -> brick 8 node el.]
%
%  INPUT
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    ed:     element displacement matrix
%
%    plotpar=[  linetype, linecolor, nodemark] 
%
%             linetype=1 -> solid    linecolor=1 -> black
%                      2 -> dashed             2 -> blue
%                      3 -> dotted             3 -> magenta
%                                              4 -> red
%             nodemark=1 -> circle       
%                      2 -> star              
%                      0 -> no mark 
%
%    magnfac:  magnification factor for displacements 
%            
%    Rem. Default is auto magnification and solid white lines with 
%         circles at nodes if magnfac and plotpar is left out
%-------------------------------------------------------------

% LAST MODIFIED: A Olsson 2004-09-03
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==4)|(nargin==5)|(nargin==6))
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
  
 d=size(ed);
 
 if ~(d(1)==a(1))
    disp('??? Check size of displacement input arguments!')
    return 
 end
 
 ned=d(2);
 
 dxmax=max(max(ex')-min(ex'));
 dymax=max(max(ey')-min(ey')); 
 dzmax=max(max(ez')-min(ez'));
 dlmax=max([dxmax;dymax;dzmax]);
 edmax=max(max(abs(ed)));
 krel=0.1;
 
 if nargin==4; 
    plotpar=[2 1 1];
    magnfac=krel*dlmax/edmax;
 elseif nargin==5;
    magnfac=krel*dlmax/edmax;
 end
 
 [s1,s2]=pltstyle(plotpar);
 
 k=magnfac;
    
% ********** Bar or Beam elements *************
% ------- Currently treated in the same way !!! ----------
    if nen==2
       if ned==6
          x=(ex+k*ed(:,[1 4]))'; 
          y=(ey+k*ed(:,[2 5]))'; 
          z=(ez+k*ed(:,[3 6]))';  
          xc=x; yc=y; zc=z;
       elseif ned==12
          x=(ex+k*ed(:,[1 7]))'; 
          y=(ey+k*ed(:,[2 8]))'; 
          z=(ez+k*ed(:,[3 9]))';  
          xc=x; yc=y; zc=z;
       end
%**********************************************************          
    else
       disp('!!!!! Sorry, this element is currently not supported!')
       return      
    end
% ************* plot commands *******************
    axis('equal')
    hold on  
    view(3)
    plot3(xc,yc,zc,s1)   
    plot3(x,y,z,s2)
    xlabel('x'); ylabel('y'); zlabel('z');
    hold off 
%--------------------------end-------------------------------- 
