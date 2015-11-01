function [sfac]=elflux2(ex,ey,es,plotpar,sfac)
%elflux2(ex,ey,es,plotpar,sfac)
%[sfac]=elflux2(ex,ey,es,plotpar)
%[sfac]=elflux2(ex,ey,es)
%-------------------------------------------------------------
% PURPOSE 
%   Display element flow arrows for a number of 2D scalar elements of 
%   the same type. To display elements and nodes, use eldraw2. 
%   Supported elements are:
%
%     1) -> triangular 3 node el.    2) -> quadrilateral 4 node el. 
%
% INPUT    
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    es:     element flow matrix
%
%    plotpar=[ arrowtype, arrowcolor]
%
%        arrowtype=1 -> solid       arrowcolor=1 -> black
%                  2 -> dashed                 2 -> blue
%                  3 -> dotted                 3 -> magenta
%                                              4 -> red
%        
%    sfac: scalc factor  =  arrowlength / element flow
%
%    Remark: If sfac and plotpar are left out the default is 
%            auto magnification and solid black arrows.
%         
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2004-10-01
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
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

 if ~(c(1)==a(1))
    disp('??? Check size of flow input argument!')
    disp('One row for each element, i.e the mean flow in x- and y-directions !') 
    return 
 end
 
 ned=c(2); 
 
% if ned~=nen ; 
%    disp('??? This function should be used for scalar problems!')
% end

 dxmax=max(max(ex')-min(ex')); 
 dymax=max(max(ey')-min(ey'));
 lm=sum(sqrt(dxmax.^2+dymax.^2))/nel;
 
 qm=sum(sqrt(sum((es').^2)))/nel;

 krel=0.8;
 
 if nargin==3; 
    plotpar=[1 1];
    sfac=lm*krel/qm;
 elseif nargin==4;
    sfac=lm*krel/qm;
 end
% *****************************************************************
 if ~((nen==3)|(nen==4))  
     disp('Sorry, this element is currently not supported!') 
     return 
 else
% ************* plot commands ******************* 
    plotprop=[plotpar 0];
    s1=pltstyle(plotprop);
 
    q=sqrt(sum((es').^2));
    la=sfac*q; fi=atan2(es(:,2),es(:,1)); loc=zeros(size(la));
  
    x0=sum(ex')/nen; y0=sum(ey')/nen;
 
    xyl0=[-0.5 0.35 0.35 0.5 0.35 0.35;
         0   0  -0.05  0  0.05  0  ];
    nar=length(x0);

    for i=1:nar
      xyl=xyl0-0.5*loc(i)*[1 1 1 1 1 1;
                           0 0 0 0 0 0];
      v=fi(i); L=la(i);

      xy=L*[ cos(v)  -sin(v) ;
             sin(v)   cos(v) ]*xyl+[x0(i)*ones(1,6);
                                    y0(i)*ones(1,6)];
      hold on
      plot( xy(1,:)',xy(2,:)',s1)
      hold off
    end;
 end        
%--------------------------end--------------------------------
