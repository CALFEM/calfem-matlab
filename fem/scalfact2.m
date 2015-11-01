function [sfac]=scalfact2(ex,ey,ed,rat)
%[sfac]=scalfact2(ex,ey,ed,rat)
%[sfac]=scalfact2(ex,ey,ed)
%-------------------------------------------------------------
% PURPOSE 
%   Determine scale factor for drawing computational results, such as 
%   displacements, section forces or flux.
%
% INPUT
%    ex,ey:  element node coordinates
%                   
%    ed:     element displacement matrix or section force matrix
%
%    rat: relation between illustrated quantity and element size. 
%    If not specified, 0.2 is used.
%    
%-------------------------------------------------------------

% LAST MODIFIED: O Dahlblom  2004-09-15
% Copyright (c)  Division of Structural Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
%
 if ~((nargin==3)|(nargin==4))
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
 
 dxmax=max(max(ex')-min(ex')); dymax=max(max(ey')-min(ey'));
 dlmax=max(dxmax,dymax);
 edmax=max(max(abs(ed)));

 if nargin==3
    k=0.2; 
 else
    k=rat;
 end
 sfac=k*dlmax/edmax;
 
%--------------------------end--------------------------------
