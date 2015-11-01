function [sfac]=eldisp2(ex,ey,ed,plotpar,sfac)
%eldisp2(ex,ey,ed,plotpar,sfac)
%[sfac]=eldisp2(ex,ey,ed,plotpar)
%[sfac]=eldisp2(ex,ey,ed)
%-------------------------------------------------------------
% PURPOSE 
%   Draw the deformed 2D mesh for a number of elements of 
%   the same type. Supported elements are:
% 
%           1) -> bar element              2) -> beam el.  
%           3) -> triangular 3 node el.    4) -> quadrilateral 4 node el. 
%           5) -> 8-node isopar. element
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
%    sfac:  scale factor for displacements 
%            
%    Rem. Default if sfac and plotpar is left out is auto magnification 
%         and dashed black lines with circles at nodes -> plotpar=[2 1 1] 
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom 2004-10-01
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==3)|(nargin==4)|(nargin==5))
    disp('??? Wrong number of input arguments!')
    return
 end
 
 a=size(ex); b=size(ey);
 
 if (a-b)==[0 0]
    nen=a(2); 
 else
    disp('??? Check size of coordinate input arguments!')
    return
 end
 
 c=size(ed);
 
 if ~(c(1)==a(1))
    disp('??? Check size of displacement input arguments!')
    return 
 end
 
 ned=c(2);
 
 dxmax=max(max(ex')-min(ex')); dymax=max(max(ey')-min(ey'));
 dlmax=max(dxmax,dymax);
 edmax=max(max(abs(ed)));
 krel=0.1;
 
 if nargin==3; 
      plotpar=[2 1 1]; sfac=krel*dlmax/edmax;
 elseif nargin==4
      sfac=krel*dlmax/edmax;
 end
 
 [s1,s2]=pltstyle(plotpar);
 k=sfac;

% ********** Bar or Beam elements *************
    if nen==2 
       if ned==4  % -----------  Bar elements -------------
          x=(ex+k*ed(:,[1 3]))'; 
          y=(ey+k*ed(:,[2 4]))';
          xc=x;
          yc=y;
       elseif ned==6  % -------- Beam elements ------------
          x=(ex+k*ed(:,[1 4]))'; 
          y=(ey+k*ed(:,[2 5]))';
          [exc,eyc]=beam2crd(ex,ey,ed,k);
          xc=exc';
          yc=eyc';
       end  
% ********* 2D triangular elements ************
    elseif nen==3
       x=(ex+k*ed(:,[1 3 5]))';   
       y=(ey+k*ed(:,[2 4 6]))';
       xc=[x; x(1,:)];
       yc=[y; y(1,:)];
       
% ********* 2D quadrilateral elements *********
    elseif nen==4 
       x=(ex+k*ed(:,[1 3 5 7]))'; 
       y=(ey+k*ed(:,[2 4 6 8]))';
       xc=[x; x(1,:)];
       yc=[y; y(1,:)];  
 % ********* 2D 8-node quadratic elements *********
    elseif nen==8 
       x=(ex+k*ed(:,[1 3 5 7 9 11 13 15])); 
       y=(ey+k*ed(:,[2 4 6 8 10 12 14 16]));
%       xc=[x(1); x(5); x(2); x(6); x(3); x(7); x(4); x(8);x(1)];
%       yc=[y(1); y(5); y(2); y(6); y(3); y(7); y(4); y(8);y(1)];     
%
% isoparametric elements
%
    t=-1;
    n=0;
    for s=-1:0.4:1
      n=n+1;
      N1=-1/4*(1-t)*(1-s)*(1+t+s);
      N2=-1/4*(1+t)*(1-s)*(1-t+s);
      N3=-1/4*(1+t)*(1+s)*(1-t-s);
      N4=-1/4*(1-t)*(1+s)*(1+t-s);
      N5=1/2*(1-t*t)*(1-s);
      N6=1/2*(1+t)*(1-s*s);
      N7=1/2*(1-t*t)*(1+s);
      N8=1/2*(1-t)*(1-s*s);
      N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

      x1(n,:)=N*x';
      y1(n,:)=N*y';
    end;
    xc=[xc x1];
    yc=[yc y1];
    clear x1
    clear y1
%
    s=1;
    n=0;
    for t=-1:0.4:1
      n=n+1;
      N1=-1/4*(1-t)*(1-s)*(1+t+s);
      N2=-1/4*(1+t)*(1-s)*(1-t+s);
      N3=-1/4*(1+t)*(1+s)*(1-t-s);
      N4=-1/4*(1-t)*(1+s)*(1+t-s);
      N5=1/2*(1-t*t)*(1-s);
      N6=1/2*(1+t)*(1-s*s);
      N7=1/2*(1-t*t)*(1+s);
      N8=1/2*(1-t)*(1-s*s);
      N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

      x1(n,:)=N*x';
      y1(n,:)=N*y';
    end;
    xc=[xc x1];
    yc=[yc y1];
    clear x1
    clear y1
%
    t=1;
    n=0;
    for s=1:-0.4:-1
      n=n+1;
      N1=-1/4*(1-t)*(1-s)*(1+t+s);
      N2=-1/4*(1+t)*(1-s)*(1-t+s);
      N3=-1/4*(1+t)*(1+s)*(1-t-s);
      N4=-1/4*(1-t)*(1+s)*(1+t-s);
      N5=1/2*(1-t*t)*(1-s);
      N6=1/2*(1+t)*(1-s*s);
      N7=1/2*(1-t*t)*(1+s);
      N8=1/2*(1-t)*(1-s*s);
      N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

      x1(n,:)=N*x';
      y1(n,:)=N*y';
    end;
    xc=[xc x1];
    yc=[yc y1];
    clear x1
    clear y1
%
    s=-1;
    n=0;
    for t=1:-0.4:-1
      n=n+1;
      N1=-1/4*(1-t)*(1-s)*(1+t+s);
      N2=-1/4*(1+t)*(1-s)*(1-t+s);
      N3=-1/4*(1+t)*(1+s)*(1-t-s);
      N4=-1/4*(1-t)*(1+s)*(1+t-s);
      N5=1/2*(1-t*t)*(1-s);
      N6=1/2*(1+t)*(1-s*s);
      N7=1/2*(1-t*t)*(1+s);
      N8=1/2*(1-t)*(1-s*s);
      N=[ N1, N2, N3 ,N4, N5, N6, N7, N8 ];

      x1(n,:)=N*x';
      y1(n,:)=N*y';
    end;
    xc=[xc x1];
    yc=[yc y1];
    clear x1
    clear y1
%
%**********************************************************       
    else
       disp('Sorry, this element is currently not supported!') 
       return
    end
% ************* plot commands *******************
    hold on
    axis equal
    plot(xc,yc,s1)
    plot(x,y,s2)
    hold off
%--------------------------end--------------------------------
