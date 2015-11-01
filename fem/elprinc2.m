function [sfac]=elprinc2(ex,ey,es,plotpar,sfac)
%elprinc2(ex,ey,es,plotpar,sfac)
%[sfac]=elprinc2(ex,ey,es,plotpar)
%[sfac]=elprinc2(ex,ey,es)
%-------------------------------------------------------------
% PURPOSE 
%   Display element principal stresses as arrows for a number of  
%   2D structural elements of the same type. To display elements and 
%   nodes, use eldraw2. Supported elements are:
%
%     1) -> triangular 3 node el.    2) -> quadrilateral 4 node el. 
%
% INPUT    
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    es:     element stress matrix
%
%    plotpar=[  arrowtype, arrowcolor]
%
%        arrowtype=1 -> solid       arrowcolor=1 -> black
%                  2 -> dashed                 2 -> blue
%                  3 -> dotted                 3 -> magenta
%                                              4 -> red
%        
%    sfac:  =  arrowlength / max element principal
%
%    Rem. Default is auto magnification and solid white arrows if sfac
%         and plotpar are left out.
%         
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
 
 a=size(ex); b=size(ey); c=size(es);
 
 if (a-b)==[0 0]
     nel=a(1);nen=a(2); 
 else
    disp('??? Check size of coordinate input arguments!') 
    return
 end
 
 if ~((nen==3)|(nen==4))  
    disp('Sorry, this element is currently not supported!') 
    return
 end    
 
 if ~(c(1)==a(1))
    disp('??? Check size of stress input argument!')
    disp('One row for each element, i.e the mean stress sig x, -y and tau xy !') 
    return 
 end
% 
% if ned~=nen ; 
%    disp('??? This function should be used for structure problems!')
% end
%
% ******** calculation of principal stresses and directions **********

 sigx=es(:,1); sigy=es(:,2); tau=es(:,3);
 
 ds=(sigx-sigy)/2; R=sqrt(ds.^2+tau.^2);
 
 sig1=(sigx+sigy)/2+R; sig2=(sigx+sigy)/2-R; alfa=atan2(tau,ds)/2;
 
% ************************************************************************
 dxmax=max(max(ex')-min(ex'));  dymax=max(max(ey')-min(ey'));
 
 lm=sum(sqrt(dxmax.^2+dymax.^2))/nel;
 
 sig1m=sum(sig1)/nel; sig2m=sum(sig2)/nel; sigm=max(sig1m,sig2m);

 krel=0.8;
 
 if nargin==3; 
    plotpar=[1 1];   sfac=lm*krel/sigm;
 elseif nargin==4;
    sfac=lm*krel/sigm;
 end

 plotpar=[plotpar 0]; s1=pltstyle(plotpar);
 
 x0=sum(ex')/nen; y0=sum(ey')/nen;
 
 la1=0.5*sfac*sig1;  la2=0.5*sfac*sig2;
 
 qrt=(pi/2).*ones(size(alfa));
 
% ************* plot commands *******************

 loc1=-sign(sig1);  loc2=-sign(sig2);

 fi=alfa;       arrow2(x0,y0,la1,fi,loc1,s1);
 
 fi=alfa+2*qrt;   arrow2(x0,y0,la1,fi,loc1,s1);      

 fi=alfa+qrt;   arrow2(x0,y0,la2,fi,loc2,s1);
 
 fi=alfa+3*qrt;   arrow2(x0,y0,la2,fi,loc2,s1);      
%--------------------------end--------------------------------
