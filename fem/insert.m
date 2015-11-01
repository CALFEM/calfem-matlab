 function [f]=insert(edof,f,fe)
% [f]=insert(edof,f,fe)
%-------------------------------------------------------------
% PURPOSE
%  Assembel fe into the global force vector f 
%  according to the topology matrix edof.
%
% INPUT: edof:  topology matrix
%        f   :  the global force vector
%        fe  :  element force vector
%
% OUTPUT:  f :  the new global force vector
%-------------------------------------------------------------

% LAST MODIFIED: M Ristinmaa 1993-10-30
% Copyright (c)  Division of Structural Mechanics and
%                Department of Solid Mechanics.
%                Lund Institute of Technology
%-------------------------------------------------------------
    [nie,n]=size(edof);
    t=edof(:,2:n);
    for i = 1:nie
         f(t(i,:))=f(t(i,:))+fe(i,:)';
    end
%--------------------------end--------------------------------
