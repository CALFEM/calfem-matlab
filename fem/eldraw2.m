function eldraw2(ex,ey,plotpar,elnum)
%eldraw2(ex,ey,plotpar,elnum)
%eldraw2(ex,ey,plotpar)
%eldraw2(ex,ey)
%-------------------------------------------------------------
% PURPOSE 
%   Draw the undeformed 2D mesh for a number of elements of 
%   the same type. Supported elements are:
%
%   1) -> bar element              2) -> beam el.  
%   3) -> triangular 3 node el.    4) -> quadrilateral 4 node el. 
%   5) -> 8-node isopar. elemen
%
% INPUT  
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements
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
%    elnum=edof(:,1) ; i.e. the first column in the topology matrix
%         
%    Rem. Default is solid white lines with circles at nodes.
%         
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom 2004-09-28
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==2)|(nargin==3)|(nargin==4))
    disp('??? Wrong number of input arguments!')
    return
 end

 a=size(ex); b=size(ey);
 
 if (a-b)==[0 0]
     nel=a(1);nen=a(2);
 else
    disp('??? Check size of coordinate input arguments!')
    return
 end
 if nargin==2; 
      plotpar=[1 1 1]; 
 end
 [s1,s2]=pltstyle(plotpar);
 
% ************************************************** 
% ************* plot coordinates *******************
% ************************************************** 
 x0=sum(ex')/nen; y0=sum(ey')/nen;
    
% ********** Bar or Beam elements *************
 if nen==2 
    x=ex' ;
    y=ey';  
    xc=x ;yc=y; 
        
% ********* 2D triangular elements ************
 elseif nen==3 
    x=ex' ;
    y=ey';
    xc=[x ; x(1,:)];  yc=[y ; y(1,:)];
          
% ********* 2D quadrilateral elements *********
 elseif nen==4 
    x=ex' ;
    y=ey';   
    xc=[x ; x(1,:)]; yc=[y ; y(1,:)];
% ********* 2D 8 node quadratic elements *********
 elseif nen==8
    x=ex;
    y=ey;    
%   xc=[x(1);x(5);x(2);x(6);x(3);x(7);x(4);x(8);x(1)];
%  yc=[y(1);y(5);y(2);y(6);y(3);y(7);y(4);y(8);y(1)];
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
%**********************************************************       
 else
    disp('!!!! Sorry, this element is currently not supported!') 
    return
 end
% ************************************************** 
% **************** plot commands *******************
% ************************************************** 
 hold on
 axis equal
 plot(xc,yc,s1) 
 plot(x,y,s2)
 if nargin==4
    for i=1:nel
        h=text(x0(i),y0(i),int2str(elnum(i)));
        set(h,'fontsize',8);
    end
 end
 %xlabel('x'); ylabel('y');
 hold off 
%--------------------------end--------------------------------
