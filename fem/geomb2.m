function geomb2(ex,ey,plotpar,elnum)
%geomb2(ex,ey,plotpar,elnum)
%geomb2(ex,ey,plotpar)
%geomb2(ex,ey)
%--------------------------------------------------------------------------
% PURPOSE 
%   Draw geometry for two dimensional bar or beam element
%
% INPUT  ex = [x1 x2]
%        ey = [y1 y2]     element node coordinates
%
%    plotpar=[linetype, linecolour, nodemark]
% 
%             linetype=1 -> solid    linecolor=1 -> black
%                      2 -> dashed             2 -> blue
%                      3 -> dotted             3 -> magenta
%                                              4 -> red
%
%             nodemark=0 -> no mark 
%                      1 -> circle       
%                      2 -> star              
%                      3 -> point 
%              
%    elnum=[element number]
%         
%    Rem. Default is solid black line with circles at nodes.
%    -> plotpar=[1 1 1]     
%--------------------------------------------------------------------------

% LAST MODIFIED: O Dahlblom 2015-11-18
% Copyright (c)  Division of Structural Mechanics and
%                Division of Solid Mechanics.
%                Lund University
%--------------------------------------------------------------------------
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
 
% ************* plot coordinates *******************

x0=sum(ex')/nen; y0=sum(ey')/nen;
    
 if nen==2 
    x=ex' ;
    y=ey';  
    xc=x ;yc=y; 
 end
 
 % **************** plot commands *******************

 hold on
 axis equal
 plot(xc,yc,s1) 
 if plotpar(3)~=0
     plot(x,y,s2)
 end    
 if nargin==4
    for i=1:nel
        h=text(x0(i),y0(i),int2str(elnum(i)));
        set(h,'fontsize',8);
    end
 end
 %xlabel('x'); ylabel('y');
%--------------------------end---------------------------------------------
  hold off 
 end