function eliso2(ex,ey,ed,isov,plotpar)
%eliso2(ex,ey,ed,isov)
%eliso2(ex,ey,ed,isov,plotpar)
%-------------------------------------------------------------
% PURPOSE 
%   Display element iso lines for a number of 2D scalar elements of 
%   the same type. To display elements and nodes, use eldraw2. 
%   Supported elements are:
%
%     1) -> triangular 3 node el.    2) -> quadrilateral 4 node el. 
%
% INPUT    
%    ex,ey:.......... nen:   number of element nodes
%                     nel:   number of elements   
%    ed:     element field variable matrix
%
%    isov:  vector containing the values of the isolines
%           if there only is one value in isov, this value is 
%           iterpreted as the number of isolines which should be
%           drawn. The distance betwen the lines will be equal.
%  
%    plotpar=[  linetype, linecolor]
%
%        linetype= 1 -> solid       linecolor= 1 -> black  
%                  2 -> dashed                 2 -> blue
%                  3 -> dotted                 3 -> magenta
%                                              4 -> red
%   
%        textfcn = 0 ->  the values of the lines will not be printed
%                  1 ->  the values of the lines will be printed at the lines     
%                  2 ->  the values of the lines will be printed where the 
%                        cursor is placed and clicked. (gtext is used)    
%
%         
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom 2004-10-01
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==4)|(nargin==5))
    disp('??? Wrong number of input arguments!')
    return
 end
 
 a=size(ex); b=size(ey); c=size(ed);niv=length(isov);
 
 if (a-b)==[0 0]
     nel=a(1);nen=a(2); 
 else
    disp('??? Check size of coordinate input arguments!') 
    return
 end

 if ~(c(1)==a(1))
    disp('??? Check size of Ed matrix!')
    disp('One row for each element, i.e the mean flow in x- and y-directions !') 
    return 
 end
 
 ned=c(2); 
 
% if ned~=nen ; 
%    disp('??? This function should be used for scalar problems!')
% end

if niv==1
   niv=isov(1);
   isomax=max(max(ed));
   isomin=min(min(ed));
   isodelta=(isomax-isomin)/(niv+1);
   for i=1:niv
      isov(i)=isomin+ i*isodelta;
   end
end
if nargin==5
 if length(plotpar)==2
    plotpar(3)=0;
 end   
end
 if nargin==4; 
    plotpar=[1 1 0];
 end
 plotpar=[plotpar(1:2) 0];
 s1=pltstyle(plotpar);
 %axis([min(min(ex)) max(max(ex)) min(min(ey)) max(max(ey))]);
% *****************************************************************
 if ~((nen==3)|(nen==4))  
     disp('Sorry, this element is currently not supported!') 
     return 
 else
 %
 ex(:,nen+1)=ex(:,1); ey(:,nen+1)=ey(:,1);ed(:,nen+1)=ed(:,1);
 test=0;
 for k=1:niv     
    for i=1:nel
        l=0;fix=0;fiy=0;
        for j=1:nen 
          if ed(i,j)<isov(k) 
             if ed(i,j+1)>=isov(k)
               test=99;
             end
          elseif ed(i,j)>isov(k)
             if ed(i,j+1)<=isov(k)
               test=99;
             end
          elseif ed(i,j)<=isov(k) 
             if ed(i,j+1)>isov(k)
               test=99;
             end
          elseif ed(i,j)>=isov(k)
             if ed(i,j+1)<isov(k)
               test=99;
             end
          end
          if test==99
            l=l+1;
            fix(l)=ex(i,j)+(ex(i,j+1)-ex(i,j))*(isov(k)-ed(i,j))....
             /(ed(i,j+1)-ed(i,j));
            fiy(l)=ey(i,j)+(ey(i,j+1)-ey(i,j))*(isov(k)-ed(i,j))....
             /(ed(i,j+1)-ed(i,j));
          end
          test=0;
        end
        if fix(1)>0.0000001
         hold on
         plot(fix,fiy,s1)
         hold off
         fixsa=fix(l);fiysa=fiy(l);
        end
    end 
    if plotpar(3)==1
       text(fixsa,fiysa,sprintf('%6.3f',isov(k)))
    elseif plotpar(3)==2
       disp('Place the cursor where you want the isovalue')
       gtext(sprintf('%6.3f',isov(k)))
    end
end
end
%--------------------------end--------------------------------
