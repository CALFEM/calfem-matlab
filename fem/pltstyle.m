function [s1,s2]=pltstyle(plotpar)
%-------------------------------------------------------------
% PURPOSE 
%   Define define linetype,linecolor and markertype character codes. 
%
% INPUT 
%    plotpar=[ linetype, linecolor, nodemark ]
% 
%             linetype=1 -> solid    linecolor=1 -> black
%                      2 -> dashed             2 -> blue
%                      3 -> dotted             3 -> magenta
%                                              4 -> red
%
%             nodemark=1 -> circle       
%                      2 -> star               
%                      0 -> no mark             
% OUTPUT
%     s1: linetype and color for mesh lines
%     s2: type and color for node markers
%-------------------------------------------------------------

% LAST MODIFIED: Ola Dahlblom 2004-09-15
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if plotpar(1)==1 ; s1='-';
 elseif plotpar(1)==2 ; s1='--';
 elseif plotpar(1)==3 ; s1=':';
 else disp('??? Error in variable plotpar(1)!');
      return;
 end
 
 if plotpar(2)==1 ; s1=[s1,'k'];
 elseif plotpar(2)==2 ; s1=[s1,'b'];
 elseif plotpar(2)==3 ; s1=[s1,'m'];
 elseif plotpar(2)==4 ; s1=[s1,'r'];
 else disp('??? Error in variable plotpar(2)!');
      return;
 end

 if plotpar(3)==1 ; s2='ko';
 elseif plotpar(3)==2 ; s2='k*';
 elseif plotpar(3)==0 ; s2='k.';
 else disp('??? Error in variable plotpar(3)!');
      return;
 end

%--------------------------end-------------------------------- 
 
